# Copyright (c) 2014-2018: Elodie Laine, Hugues Richard, Adel Ait-hamlat
# and Diego Javier Zea.
# This code is part of the phylosofs package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

# Modele d'inference de phylogenies de transcrits
# Pour etudier l'apparition et la fixation
# d'evenements d'epissage au cours de l'evolution

# all necessary functions

from __future__ import division
import time
import random
import subprocess
import copy
import os
import networkx as nx
import numpy as np
import math as M
import pickle as pk
import drawForestTranscripts as df
import warnings
import nx_utils
##################################################

#            FONCTIONS INPUT/OUTPUT              #

##################################################


# fonction prend en parametre un object python o et un nom de fichier s
# et va creer ou ecraser le fichier de nom s avec l'objet ecrit en binaire
def pickPrint(o, s):
    with open(s, "wb") as f:
        pk.dump(o, f)

# prochaine fonction inverse de la precedente qui cree un objet python
# a partir du fichier de nom s ou il a ete ecrit en binaire


def pickRead(s):
    with open(s, "rb") as f:
        return (pk.load(f))

# print all donfigurations as  forests of trees
# each transcript family (tree) is indicated by a unique color
# a file with extension .info is also created
# where each ancestral and current transctipt is labelled


def drawForest(filename, outputDir):
    with open(os.path.join(outputDir, filename+".pk"), 'rb') as f:
        transU = pk.load(f)
        count = 0
        for i in transU:
            df.ForestToDot(i, os.path.join(outputDir, filename), count)
            subprocess.call(["dot", "-Tpdf", "-o",
                             os.path.join(outputDir, filename + "_config" +
                                          str(count) + ".pdf"),
                             os.path.join(outputDir, filename + "_config" +
                                          str(count) + ".dot")])
            count = count + 1

# write the output results  the text files
# and print all configurations as forests of trees


def writeOutput(res, exSt, SUFF, costs, AllExons, outputDir):
    print "\n>>> writing output..."
    tree = res[0]
    res2 = tree_cost(tree, 10000, costs)
    confs = res2[1]
    treeCost = res2[2]
    Res = []
    pickPrint(tree, os.path.join(outputDir, "tree" + SUFF + ".pk"))
    Tree = transTree(tree, AllExons)
    with open(os.path.join(outputDir, "final" + SUFF + ".txt"), 'w') as ff:
        ff.write(
            "\n~~~~~~~~   Phylosofs results  :\n\n"
            "   ** exon state for ancestral nodes : \n\n")
        ff.write(nx_utils.str_nodes(exSt))
        ff.write("\n\n   ** Best topology found : \n\n")
        ff.write(nx_utils.str_nodes(tree))
        ff.write("\n\n   ** The transcript tree for this topology : \n\n")
        ff.write(nx_utils.str_nodes(Tree))
        ff.write("\n\n   ** Score for this topology = " +
                 str(treeCost) + "\n\n")
        ff.write("   ** Execution time : " + str(res[2]) +
                 " ; number of cuts : " + str(res[1]) + "\n")
    for c in confs:
        Res.append(elagTree(tree, c, AllExons))
    pickPrint(Res, os.path.join(outputDir, "solution" + SUFF + ".pk"))
    drawForest("solution" + SUFF, outputDir)

    ##################################################

    #              FONCTIONS GENERALES               #

    ##################################################


# Un transcrit est une chaine de caractere
# Chez une espece, a un gene correspond une liste de transcrits

# Notes on the data structure used to represent transcripts' phylogenies
# (forests of phylogenetic trees):
# - DiGraph from NetworkX package
# - each node represents an ancestral (internal nodes) or
# current(leaves) gene of a species
# - a number of possible and equivalent configurations
# is given for each node(attribute of the node: configuration)
# there is a dependency between the numbers of configurations
# of the different nodes(the attributes rConf and lConf
# of the considered configuration of the node give
# the indices of the associated configurations
# in the left and right children, the attribute parConf of
# inherited gives the associated configuration of the parent)
# - each node is divided into as many subnodes as transcripts
# which numbers are given by the attributes ln, rn and bn
# - the exon compositions of the transcripts are given
# by the attributes rft, Bt and lft, each exon is represented
# by a list of costs corresponding to [0/0 or 1/0, 1/1, 2/1]
# - the assignments are stored in the attribute inherited
# parConf correspond to the configurations of the parent(s)
# we go up to the point where we find an ancestral binary subnode
# the variables binary, right and left gives the assigment for the
# binary, right and left subnodes of the considered node(the child)
# the second index in each tuple corresponds to the index of the transcript
# in Bt, rft and lft, the first index corresponds to the index of the
# transcript in Bt of the retained ancestral node

# La fonction exonState prend une liste de transcrits
# et renvoie l'etat de chque exon dans cet ensmble

# exonState : string' list -> {0,1,2}' list

def exonState(g, AllExons):
    res = []
    g_len = len(g)
    for e in AllExons:
        tmp = 0
        for t in g:
            tmp = tmp + (e in t)
        if tmp == 0:
            res.append(0)
        elif tmp == g_len:
            res.append(2)
        else:
            res.append(1)
    return res


# les noeuds d'un arbre sont des entier. Leaves prend
# un arbre oriente et renvoie la liste des feuilles de cet arbre

# leaves : nx.DiGraph() -> int list

def leaves(t):
    res = set([])
    for n in t.nodes():
        if nx_utils.successors(t, n) == []:
            res.add(n)
    return res


# parents est utile pour toute fonction qui remonte recursivement dans un arbre
# elle renvoie a partir d'un ensemble de noeuds l'ensemble de noeuds au
# niveau superieur cette fonction prend en argument :
#   -gt  : un arbre oriente
#   -ens : une liste de deux ensembles :
#   - ens[0] : la liste sur la quelle on itere pour les remontees classiques
#   - ens[1] : noeuds dont tous les enfants sont dans ens[0]


# parents : nx.DiGraph , [int' set , int' set] ->  [int' set , int' set]

def parents(gt, ens):
    res1 = set([])
    res2 = set([])
    for i in ens[1]:
        if nx_utils.predecessors(gt, i) != []:
            try:
                tmp0 = nx_utils.predecessors(gt, i)[0]
                tmp = {tmp0}
            except IndexError:
                tmp = set([])
            if len(nx_utils.successors(gt, tmp0)) == 2:
                l, r = nx_utils.successors(gt, tmp0)
                if (l in ens[1] and r in ens[1]):
                    res1 = res1 | tmp
                    res2 = res2 | tmp
                else:
                    res2 = res2 | {i}
            else:
                res1 = res1 | tmp
                res2 = res2 | tmp
    return([res1, res2])


# STRUCTURES DE DONNEES :

# INPUTS :

# Arbre phylogenetique des especes etudiees pour le gene donne
# (chaque noeud de cet arbre est une feuille ou a deux descendants.

# Les transcrits observes pour ce gene chez chaque espece

# Un arbre de genes est un arbre oriente dont les fuilles ont un label 'trans'
# contenant la liste des transcrtis observes chez l'espece correspondante

# Une foret de topologies associee a un arbre de gene est une arbre oriente
# Cahque noeud interne de cet arbre contient les labels suivants :
# bn (Binary (sub)Nodes) :
# le nombre de noeuds binaires. Un noeud binaire est un transcrit de l'espece
# ancestral correspondante qui a deux transcrits descendants,
# un dans chacune des deux especes filles.
# ln (Left(sub)Nodes) :
# le nombre de noeuds gauche, i.e un transcrit ayant un seul descendant
#
# Un arbre de transcrits est arbre oriente construit a partir
# d'un arbre de gene
# Chaque noeud d'un arbre de transcrits contient :
#

# aleatopo prend un arbre de gene gt et renvoie une
# foret aleatoire de topologies coherente avec gt

# aleatopo : nx.DiGraph -> nx.DiGraph

def aleaTopo(gt, minBNTree):

    tt = nx.DiGraph()
    tt.add_edges_from(gt.edges())
    tmp = leaves(gt)
    ens = [tmp, tmp]

    for i in tmp:
        tt.node[i]['trans'] = gt.node[i]['trans']
        tt.node[i]['n'] = len(gt.node[i]['trans'])

    ens = parents(gt, ens)

    # l'algo part d'en bas
    while ens[0] != set([]):
        # print ens[0]
        for e in ens[0]:
            l, r = nx_utils.successors(gt, e)
            # print "je suis "+str(e)+" "+str(l)+" "+str(r)
            ml = tt.node[l]['n']
            mr = tt.node[r]['n']
            # print ml,mr
            # print tt.node[l]['trans']
            # print tt.node[r]['trans']
            minBN = minBNTree.node[e]['minBN']
            # print minBNTree.nodes()
            # print type(minBN)
            # print minBN
            bn = random.randint(minBN, min(ml, mr))
            if nx_utils.predecessors(gt, e) == []:
                tt.node[e]['bn'] = bn
                tt.node[e]['n'] = bn
                tt.node[e]['ln'] = 0
                tt.node[e]['rn'] = 0
            else:
                tt.node[e]['bn'] = bn
                ln = random.randint(bn == 0, ml-bn)
                tt.node[e]['ln'] = ln
                rn = random.randint(bn == 0, mr-bn)
                tt.node[e]['rn'] = rn
                tt.node[e]['n'] = bn+ln+rn

        ens = parents(gt, ens)
    return (tt)


# next function takes a list l
# and returns the sum of the n mins of l
def minScore(n, l):
    if n <= 0:
        return(0)
    else:
        tmp = min(l)
        indTmp = l.index(tmp)
        l.pop(indTmp)
        return (tmp + minScore(n-1, l))


# fonction suivante renvoie l'indice du min de l
def argMin(l):
    if l != []:
        res = 0
        tmp = l[0]
        for i in range(1, len(l)):
            if l[i] < tmp:
                res = i
        return(res)
    else:
        print("problem : argmin sur list vide\n")
        return(-1)

# fonction suivante renvoie l'indice ou les indices du min de l


def argMinMul(l):
    if l != []:
        res = []
        m = min(l)
        for i in range(len(l)):
            if l[i] == m:
                res.append(i)
        return(res)
    else:
        print("problem : argmin sur list vide\n")
        return(-1)


# pour l'instant on supposera que la racine de tout arbre = 1
def racine(t):
    return(1)

# fonction qui renvoie tous les indices de l[ind:] ou n est present
#         indexes : elmt' list , int , elmt -> int' list


def indices(l, ind, n):
    if (l == []):
        return ([])
    else:
        try:
            i = l[ind:].index(n)+ind
            return ([i] + indices(l, i+1, n))
        except ValueError:
            return([])

            ##################################################

            #       ETAT DES EXONS AUX NOEUDS INTERNES       #

            ##################################################


# sank is a fun which takes a tree t and a node n
# and will generate the cost as sankoff algo of the node n.
# children of n in t must have a key Cost as np.array
# exS is a list of booleans of length of AllExons
# for each exon e , exS[e]==1 iff e has to be present
# in the transcripts of node n in repspect to
# the Dollo's parcimony.
# this is for step 1 : determine the states of the
# exons at the genome level

def sank(t, n, exS, costMat, AllExons):

    C0 = costMat[0]
    # C1_0 = costMat[1]
    # C1_1 = costMat[2]
    # C1_2 = costMat[3]
    nExons = range(len(AllExons))

    l, r = nx_utils.successors(t, n)
    res = []
    # for each exon
    for e in nExons:
        # get the costs associated to exon e in the left and right
        # children of n (l and r)
        cl = t.node[l]['Cost'][e]
        cr = t.node[r]['Cost'][e]
        tmpRes = []
        for i in [0, 1, 2]:
            tmpCl = C0[i][0] + cl[0, 0]
            tmpIl = 0
            tmpCr = C0[i][0] + cr[0, 0]
            tmpIr = 0
            for j in [1, 2]:
                if (C0[i][j] + cl[j, 0]) < tmpCl:
                    tmpCl = C0[i][j] + cl[j, 0]
                    tmpIl = j
                if (C0[i][j] + cr[j, 0]) < tmpCr:
                    tmpCr = C0[i][j] + cr[j, 0]
                    tmpIr = j
            tmpRes.append([tmpCl+tmpCr, tmpIl, tmpIr])
        if (exS[e] == 1):
            tmpRes[0] = [10000, -2, -2]
        res.append(tmpRes)
    t.node[n]['Cost'] = np.array(res)

# propagateOnes utile pour dollo


def propagateOnes(t, n, e):
    tmp = t.node[n]['dollo'][e]
    if nx_utils.successors(t, n) != []:
        l, r = nx_utils.successors(t, n)
        if tmp == 2:
            t.node[n]['dollo'][e] = 1
            propagateOnes(t, l, e)
            propagateOnes(t, r, e)

# dollo sert a forcer l'etat des exons a respecter
# la condition de dollo


def dollo(t, AllExons):
    nExons = range(len(AllExons))
    f = leaves(t)
    res = nx.DiGraph()
    res.add_nodes_from(t.nodes())
    res.add_edges_from(t.edges())
    for i in f:
        exS = exonState(t.node[i]['trans'], AllExons)
        for j in nExons:
            if exS[j] > 0:
                exS[j] = 1
        res.node[i]['dollo'] = exS
    ens = parents(t, [f, f])
    while(ens[0] != set([])):
        for e in ens[0]:
            tmp = []
            l, r = nx_utils.successors(t, e)
            exSL = res.node[l]['dollo']
            exSR = res.node[r]['dollo']
            for ee in nExons:
                ll, rr = exSL[ee], exSR[ee]
                if (ll, rr) == (0, 0):
                    tmp.append(0)
                elif ll == 0 or rr == 0:
                    tmp.append(2)
                elif ll == 1 and rr == 1:
                    tmp.append(1)
                else:
                    tmp.append(1)
                    propagateOnes(res, l, ee)
                    propagateOnes(res, r, ee)
                res.node[e]['dollo'] = tmp
        ens = parents(t, ens)
    return (res)


# sankBack takes a cost tree ct and the exon state tree st
# l the list of the possible states of node n and will change recursively
# st.node[n][extSate] for the exon e then will do the same for children

def sankBack(ct, st, n, el, AllExons):
    nExons = range(len(AllExons))
    if nx_utils.successors(ct, n) != []:
        l, r = nx_utils.successors(st, n)
        st.node[n]['exState'] = el
        tmpL = []
        tmpR = []
        for i in nExons:
            tmpL0 = []
            tmpR0 = []
            for j in el[i]:
                ll, rr = ct.node[n]['Cost'][i, j][1:]
                tmpL0.append(ll)
                tmpR0.append(rr)
            tmpL.append(tmpL0)
            tmpR.append(tmpR0)
        st.node[l]['exState'] = []
        st.node[r]['exState'] = []
        sankBack(ct, st, l, tmpL, AllExons)
        sankBack(ct, st, r, tmpR, AllExons)


# exState prend un gene tree et renvoie un esxonState tree
# ou chaque noeud contient ExState l'etat de chaque exon dans l'ensemble
# de ses transcrits

def exState(t, costMat, AllExons):
    nExons = range(len(AllExons))
    tmp = nx.DiGraph()
    # nodes = t.nodes()
    tmp.add_edges_from(t.edges())
    f = leaves(t)
    res = nx.DiGraph()
    res.add_edges_from(t.edges())
    for i in f:
        exState_of_i = exonState(t.node[i]['trans'], AllExons)
        res.node[i]['est'] = exState_of_i
        cost = []
        for e in nExons:
            costTmp = [[10000, -1, -1], [10000, -1, -1], [10000, -1, -1]]
            costTmp[exState_of_i[e]][0] = 0
            cost.append(costTmp)
        tmp.node[i]['Cost'] = np.array(cost)

    dolloMask = dollo(t, AllExons)
    ens = parents(t, [f, f])
    while(ens[0] != set([])):
        for e in ens[0]:
            exS = dolloMask.node[e]['dollo']
            sank(tmp, e, exS, costMat, AllExons)
        ens = parents(t, ens)

    r = racine(t)
    res.node[r]['exState'] = []
    rootCost = tmp.node[r]['Cost']
    minCost = rootCost.min(axis=1)[:, 0]
    tmp2 = []
    for e in nExons:
        minCost_of_e = minCost[e]
        ll = rootCost[e, :, 0].tolist()
        indMins = indices(ll, 0, minCost_of_e)
        tmp2.append(indMins)

    sankBack(tmp, res, r, tmp2, AllExons)

    return (res)

# ex_state permet de ne prendre qu'une configuration parmi toutes
# celles possibles des etats possibles des exons pour chaque noeud.
# t est le resultat de exState appliquee a un arbre de genes
# la priorite est donnee comme variable globale

# ex_state : nx.DiGraph() -> nx.DiGraph()


def ex_state(t, priority):
    tt = nx.DiGraph()
    tt.add_edges_from(t.edges())
    for i in tt.nodes():
        if nx_utils.successors(tt, i) != []:
            tt.node[i]['est'] = []
            for ii in t.node[i]['exState']:
                if priority[0] in ii:
                    tmp = priority[0]
                elif priority[1] in ii:
                    tmp = priority[1]
                else:
                    tmp = priority[2]
                tt.node[i]['est'].append(tmp)
        else:
            tt.node[i]['est'] = t.node[i]['est']
    return(tt)

    ##################################################

    #                Leaf Assignment                 #

    ##################################################


############
# Scores : #
############

def toScore(trans, exS, AllExons):
    nExons = range(len(AllExons))
    res = []
    # print "TRANS toScore --> " + str(trans)
    # print "AllExons toScore --> " + str(AllExons)
    for i in trans:
        tmp = []
        for e in nExons:
            tmpp = [10000, 10000, 10000]
            if exS[e] == 1:
                if AllExons[e] in i:
                    tmpp[1] = 0
                else:
                    tmpp[0] = 0
            else:
                tmpp[exS[e]] = 0
            tmp.append(tmpp)
        res.append(tmp)
    return (res)


# transScore renvoie le transcrit parent de
# t1 et t2 en respect avec l'etat des exons exS


def transScore(t1, t2, exS, exS_g, exS_d, costMat, AllExons):

    # C0 = costMat[0]
    C1_0 = costMat[1]
    C1_1 = costMat[2]
    C1_2 = costMat[3]
    nExons = range(len(AllExons))

    res = []
    for e in nExons:
        cl = t1[e]
        cr = t2[e]

        if exS[e] == 1:

            tmpRes = []
            # if the child exon status is different from 1 (0 or 2),
            # then take C1_1[1] (second cell) ok
            CG = C1_1[exS_g[e] != 1]
            CD = C1_1[exS_d[e] != 1]

            for i in [0, 1]:

                tmpCr = min([CG[i][k] + cl[k] for k in [0, 1, 2]])
                tmpCl = min([CD[i][k] + cr[k] for k in [0, 1, 2]])

                tmpRes.append(tmpCl+tmpCr)
            tmpRes.append(10000)

        elif exS[e] == 0:

            tmpCl = min([C1_0[i] + cl[i] for i in [0, 1, 2]])
            tmpCr = min([C1_0[i] + cr[i] for i in [0, 1, 2]])

            tmpRes = [tmpCl + tmpCr, 10000, 10000]

        else:
            CG = C1_2[exS_g[e] != 0]
            CD = C1_2[exS_d[e] != 0]

            tmpCl = min([CG[i] + cl[i] for i in [0, 1, 2]])
            tmpCr = min([CD[i] + cr[i] for i in [0, 1, 2]])

            tmpRes = [10000, 10000, tmpCl + tmpCr]
        res.append(tmpRes)

    return (res)


# le cout d'un transcrit t = somme des mins des lignes de t
def scoreTrans(t):
    res = 0
    for i in t:
        res = res + (min(i))
    return(res)


# la fonction suivante permet de mesurer la distance entre deux transcrits

def transDist(t1, t2, exS, exs_g, exs_d, costMat, AllExons):
    res = scoreTrans(transScore(t1, t2, exS, exs_g, exs_d, costMat, AllExons))
    return(res)


##############
# best affectation
#########################

# findMins sert a trouver les n (=bn) min de a
# dim est la dimension de a et forvals les valeurs interdites
# la fonction renvoie aussi argmin les indices du min de a

def findMins(n, a, dim, forVals):
    # print "finding n min of a"
    nrows, ncols = dim
    res = []
    nn = n
    countJ = 0
    countI = 0
    while [countI, countJ] in forVals:
        nn = nn + 1
        if countJ < ncols - 1:
            countJ = countJ + 1
        else:
            countI = countI + 1
            countJ = 0
    maxMin = a[countI, countJ]
    res.append(maxMin)
    argMin = (countI, countJ)
    Min = maxMin
    argMaxMin = 0
    if countJ < ncols - 1:
        countJ = countJ + 1
    else:
        countI = countI + 1
        countJ = 0

    while countI * ncols + countJ < nn:
        if [countI, countJ] not in forVals:
            tmp = a[countI, countJ]
            if tmp > maxMin:
                maxMin = tmp
            if tmp < Min:

                argMin = (countI, countJ)
                Min = tmp
            res.append(tmp)
        else:
            nn = nn + 1
        if countJ < ncols - 1:
            countJ = countJ + 1
        else:
            countI = countI + 1
            countJ = 0
    argMaxMin = res.index(max(res))
    N = ncols * nrows
    while countI * ncols + countJ < N:
        if [countI, countJ] not in forVals:
            tmp = a[countI, countJ]
            if tmp < maxMin:
                res[argMaxMin] = tmp
                maxMin = max(res)
                argMaxMin = res.index(maxMin)
            if tmp < Min:
                argMin = (countI, countJ)
                Min = tmp
        if countJ < ncols - 1:
            countJ = countJ + 1
        else:
            countI = countI + 1
            countJ = 0
    # print "found"
    return (res, argMin)


# adaptX met a jour x en fonction de
# misLines, les lignes ou colonnes deja supprimees
def adaptX(x, misLines):
    countX = x
    n = len(misLines)
    i = 0
    crossed = [False] * n
    while i < n:
        if countX >= misLines[i][0] and not crossed[i]:
            countX = countX + misLines[i][1]
            crossed[i] = True
            i = 0
        else:
            i = i + 1
    return (countX)


# adaptX2 met a jour x en fonction de
# misLines, les lignes ou colonnes deja supprimees
# ici, a la difference de plus haut, misLines ne
# contient que les indices des lignes supprimees
# nx est le nombre de lignes supprimees a partir de x
def adaptX2(x, nx, misLines):
    countX = x
    n = len(misLines)
    i = 0
    crossed = [False] * n
    while i < n:
        if countX >= misLines[i] and not crossed[i]:
            countX = countX + 1
            crossed[i] = True
            i = 0
        else:
            i = i + 1
    return (range(countX, countX+nx))


# x premiere ligne ou colonne a supprimer et
# y le nombre de lignes supprimees a partir de x

def adaptMisLines(misLines, (x, y)):
    xx = adaptX(x, misLines)
    yy = y
    TMP = []
    misLines.sort()
    for i in range(len(misLines)):
        if misLines[i][0] > xx and misLines[i][0] < xx+yy:
            TMP.append(i)
            yy = yy+misLines[i][1]

    count = 0
    for k in TMP:
        misLines.pop(k-count)
        count = count + 1

    misLines.append((xx, yy))

# inheritedFrom sert pour celle d'apres
# elle repere dans t.node[n]['configurations'][conf][inherited]
# le dictionnaire qui a pour parConf parconf et si elle le trouve
# met a jour res en concatenant les indices de Bt herites a res[0] ..


def inheritedFrom(t, n, conf, parConf, res):
    # print "iheritedFrom"
    # print t, n, conf, parConf, res
    tmp = t.node[n]['configurations'][conf]['inherited']
    i = 0
    found = False
    while i < len(tmp) and not found:
        if tmp[i]['parConf'] == parConf:
            found = True
            for k in tmp[i]['binary']:
                res[0].append(k[1])
            for k in tmp[i]['left']:
                res[1].append(k[1])
            for k in tmp[i]['right']:
                res[2].append(k[1])
        i = i + 1

# inherited prend un arbre issu de leafAssign et
# va pour un noeud n et une configuration conf et liste des
# configurations des parents parConf, la liste des listes
# d'indices de respectivemtn bn, lft, rft,  qui sont herites.


def inherited(t, n, conf, parConf):
    lst = parConf[:]
    res = [[], [], []]
    while lst != []:
        inheritedFrom(t, n, conf, lst, res)
        lst.pop(0)
    return(res)


def testParaT(n):
    res = nx.DiGraph()
    res.add_node(n)
    res.node[n]['bn'] = 3
    res.node[n]['fbn'] = 3
    res.node[n]['ln'] = 0
    res.node[n]['rn'] = 0
    res.node[n]['nlft'] = 0
    res.node[n]['nrft'] = 0
    res.node[n]['nlt'] = 0
    res.node[n]['nrt'] = 0
    res.node[n]['lftInd'] = []
    res.node[n]['rftInd'] = []
    return(res)


def genereParaTRec0(t, n, conf, parConf, res):
    tmp = inherited(t, n, conf, parConf)

    if nx_utils.successors(t, n) == []:
        res.node[n]['bn'] = t.node[n]['bn']
        res.node[n]['fbn'] = t.node[n]['bn'] - len(tmp[0])
        res.node[n]['ln'] = 0
        res.node[n]['rn'] = 0
        res.node[n]['nlft'] = 0
        res.node[n]['nrft'] = 0
        res.node[n]['nlt'] = 0
        res.node[n]['nrt'] = 0
        res.node[n]['lftInd'] = []
        res.node[n]['rftInd'] = []
    else:

        lConf = t.node[n]['configurations'][conf]['lConf']
        rConf = t.node[n]['configurations'][conf]['rConf']

        l, r = nx_utils.successors(t, n)
        res.add_edge(n, l)
        res.add_edge(n, r)

        res.node[n]['bn'] = t.node[n]['bn']
        res.node[n]['fbn'] = t.node[n]['bn'] - len(tmp[0])

        nlt = len(t.node[n]['configurations'][conf]['lftInd'])
        tmp2 = t.node[n]['ln'] - len(tmp[1])
        res.node[n]['ln'] = tmp2
        res.node[n]['nlt'] = nlt

        nrt = len(t.node[n]['configurations'][conf]['rftInd'])
        tmp3 = t.node[n]['rn'] - len(tmp[2])
        res.node[n]['rn'] = tmp3
        res.node[n]['nrt'] = nrt

        if tmp2 == 0:
            res.node[n]['lftInd'] = []
            res.node[n]['nlft'] = 0
        else:
            res.node[n]['lftInd'] = t.node[n]['configurations'][conf]['lftInd']
            res.node[n]['nlft'] = nlt - len(tmp[1])

        if tmp3 == 0:
            res.node[n]['rftInd'] = []
            res.node[n]['nrft'] = 0
        else:
            res.node[n]['rftInd'] = t.node[n]['configurations'][conf]['rftInd']
            res.node[n]['nrft'] = nrt - len(tmp[2])

        pc1 = parConf[:]
        pc1.append(conf)
        genereParaTRec(t, l, lConf, pc1, res)
        pc2 = parConf[:]
        pc2.append(conf)
        genereParaTRec(t, r, rConf, pc2, res)


# remonteTPRec utile pour remonteTP sert a basiier nlft ou nrft
# dans le cas d'un ln ou rn = 0 rencontrE en construction de paraT
# t arbre de transcrits, pt paraT, n noeud, par parent, grdpar
def remonteTPRec(t, pt, n, par, grdPar, conf, PC0, PC1, un, nft, inh):
    enf = nx_utils.successors(pt, grdPar)
    if par == enf[0]:
        if pt.node[grdPar]['ln'] > 0:
            tmp = PC0.pop()
            PC1 = [tmp] + PC1
            k = inherited(t, n, conf, PC1)
            if len(k[inh]) < un:
                pt.node[grdPar]['nlft'] = pt.node[grdPar]['nlft'] - nft
                remonteTPRec(t, pt, n, grdPar, nx_utils.predecessors(
                    pt, grdPar)[0], conf, PC0, PC1, un, nft, inh)

    else:

        if pt.node[grdPar]['rn'] > 0:
            tmp = PC0.pop()
            PC1 = [tmp] + PC1
            k = inherited(t, n, conf, PC1)
            if len(k[inh]) < un:
                pt.node[grdPar]['nrft'] = pt.node[grdPar]['nrft'] - nft
                remonteTPRec(t, pt, n, grdPar, nx_utils.predecessors(
                    pt, grdPar)[0], conf, PC0, PC1, un, nft, inh)


def remonteTP(t, pt, n, conf, parConf, un, nft, inh):

    if nx_utils.predecessors(pt, n) != []:
        par = nx_utils.predecessors(pt, n)[0]
        PC0 = parConf[:]
        PC1 = []
        remonteTPRec(t, pt, n, n, par, conf, PC0, PC1, un, nft, inh)


def genereParaTRec(t, n, conf, parConf, res):

    tmp = inherited(t, n, conf, parConf)

    if nx_utils.successors(t, n) == []:
        res.node[n]['bn'] = t.node[n]['bn']
        res.node[n]['fbn'] = t.node[n]['bn'] - len(tmp[0])
        res.node[n]['ln'] = 0
        res.node[n]['rn'] = 0
        res.node[n]['nlft'] = 0
        res.node[n]['nrft'] = 0
        res.node[n]['nlt'] = 0
        res.node[n]['nrt'] = 0
        res.node[n]['lftInd'] = []
        res.node[n]['rftInd'] = []

    else:

        lConf = t.node[n]['configurations'][conf]['lConf']
        rConf = t.node[n]['configurations'][conf]['rConf']

        l, r = nx_utils.successors(t, n)
        res.add_edge(n, l)
        res.add_edge(n, r)

        bn = t.node[n]['bn']
        res.node[n]['bn'] = bn
        res.node[n]['fbn'] = t.node[n]['bn'] - len(tmp[0])

        nlt = len(t.node[n]['configurations'][conf]['lftInd'])
        ln = t.node[n]['ln']
        tmp2 = ln - len(tmp[1])
        res.node[n]['ln'] = tmp2
        res.node[n]['nlt'] = nlt

        nrt = len(t.node[n]['configurations'][conf]['rftInd'])
        rn = t.node[n]['rn']
        tmp3 = rn - len(tmp[2])
        res.node[n]['rn'] = tmp3
        res.node[n]['nrt'] = nrt
        if tmp2 == 0:
            res.node[n]['lftInd'] = []
            res.node[n]['nlft'] = 0
            if ln > 0 and parConf != []:

                par = nx_utils.predecessors(t, n)[0]
                enf = nx_utils.successors(t, par)
                if n == enf[0]:
                    side = 'lftInd'
                    ind_I = 0
                else:
                    side = 'rftInd'
                    ind_I = 1

                tmpp = 0

                # ce qui suit craint la mort, c est a ameliorer
                for i in range(nlt):
                    if ((i + bn) in
                            t.node[par]['configurations'][parConf[-1]][side]):
                        tmpp = tmpp + 1
                    elif ((i + bn) in
                          [y[ind_I] for y in
                           t.node[par]['configurations'][parConf[-1]]['Lbn']]):
                        tmpp = tmpp + 1
                tmpp = tmpp - len(tmp[1])
                remonteTP(t, res, n, conf, parConf, ln, tmpp, 1)
        else:
            res.node[n]['lftInd'] = t.node[n]['configurations'][conf]['lftInd']
            res.node[n]['nlft'] = nlt - len(tmp[1])

        if tmp3 == 0:
            res.node[n]['rftInd'] = []
            res.node[n]['nrft'] = 0
            if rn > 0 and parConf != []:
                par = nx_utils.predecessors(t, n)[0]
                enf = nx_utils.successors(t, par)
                if n == enf[0]:
                    side = 'lftInd'
                    ind_I = 0
                else:
                    side = 'rftInd'
                    ind_I = 1

                tmpp = 0
                # ce qui suit craint la mort, c est a ameliorer
                for i in range(nrt):
                    if ((i + bn + nlt) in
                            t.node[par]['configurations'][parConf[-1]][side]):
                        tmpp = tmpp + 1
                    if (i + bn + nlt) in \
                        [y[ind_I] for y in
                            t.node[par]['configurations'][parConf[-1]]['Lbn']]:
                        tmpp = tmpp + 1

                tmpp = tmpp - len(tmp[2])
                remonteTP(t, res, n, conf, parConf, rn, tmpp, 2)
        else:
            res.node[n]['rftInd'] = t.node[n]['configurations'][conf]['rftInd']
            res.node[n]['nrft'] = nrt - len(tmp[2])

        pc1 = parConf[:]
        pc1.append(conf)
        genereParaTRec(t, l, lConf, pc1, res)
        pc2 = parConf[:]
        pc2.append(conf)
        genereParaTRec(t, r, rConf, pc2, res)


def generateParaT(t, n, conf, parConf):
    res = nx.DiGraph()
    res.add_node(n)
    genereParaTRec(t, n, conf, parConf, res)
    return(res)


def generateParaT0(t, n, conf, parConf):
    res = nx.DiGraph()
    res.add_node(n)
    genereParaTRec0(t, n, conf, parConf, res)
    return(res)


# adaptFt used for adaptForVals to send from a node n
# to the root of paraT the number of deleted lines and will
# adapt nlft/rflt of the nodes in the path
def adaptFt(n, paraT, nft):
    if nx_utils.predecessors(paraT, n) != []:
        par = nx_utils.predecessors(paraT, n)[0]
        enf = nx_utils.successors(paraT, par)
        if n == enf[0]:
            paraT.node[par]['nlft'] = paraT.node[par]['nlft'] - nft
            adaptFt(par, paraT, nft)
        else:
            paraT.node[par]['nrft'] = paraT.node[par]['nrft'] - nft
            adaptFt(par, paraT, nft)


# table of score t, x the index of line from ind0
# xt corresponding transcript in the corresponding node
# res the transcript tree beeing constructed

def adaptDistTab(t, x, xt, forVals, ind0, paraT, n, ax):

    nRows, nCols = t.shape
    bn = paraT.node[n]['bn']
    fbn = paraT.node[n]['fbn']
    ln = paraT.node[n]['ln']
    rn = paraT.node[n]['rn']
    nlft = paraT.node[n]['nlft']
    nrft = paraT.node[n]['nrft']
    nlt = paraT.node[n]['nlt']
    nrt = paraT.node[n]['nrt']
    rftInd = paraT.node[n]['rftInd']
    lftInd = paraT.node[n]['lftInd']

    # print "nRows:", nRows, "nCols:", nCols, "xt:", xt, "bn:", bn
    # print "nlt:", nlt, "nrt:", nrt

    # in case of a vector instead of an array
    if nRows == 0 or nCols == 0:

        if nx_utils.successors(paraT, n) != []:
            l, r = nx_utils.successors(paraT, n)

        if xt >= bn:  # On ne rentre jamais dans ce cas si n feuille

            tmpXT = nlt + bn

            if xt >= tmpXT:

                tmpX = nlft + fbn

                xx = x - tmpX
                xxt = rftInd[xt - tmpXT]

                paraT.node[n]['rn'] = paraT.node[n]['rn'] - 1

                tmpp = ind0 + tmpX

                if paraT.node[n]['rn'] == 0:
                    paraT.node[n]['nrft'] = 0
                    adaptFt(n, paraT, nrft)

                    return(np.array([]), tmpp, nrft)

                else:

                    return(adaptDistTab(t, xx, xxt, forVals,
                                        tmpp, paraT, r, ax))

            else:

                xx = x - fbn
                xxt = lftInd[xt - bn]

                paraT.node[n]['ln'] = paraT.node[n]['ln']-1
                tmpp = ind0 + fbn

                if paraT.node[n]['ln'] == 0:
                    paraT.node[n]['nlft'] = 0
                    adaptFt(n, paraT, nlft)

                    return (np.array([]), tmpp, nlft)

                else:

                    return(adaptDistTab(t, xx, xxt, forVals,
                                        tmpp, paraT, l, ax))

        else:

            tmp = ind0 + x
            paraT.node[n]['fbn'] = paraT.node[n]['fbn']-1

            return (np.array([]), tmp, 1)

    # in case of a real array
    else:

        if nx_utils.successors(paraT, n) != []:
            l, r = nx_utils.successors(paraT, n)

        if xt >= bn:  # On ne rentre jamais dans ce cas si n feuille

            tmpXT = nlt + bn

            if xt >= tmpXT:

                tmp = nlft + fbn
                xx = x - tmp
                xxt = rftInd[xt - nlt - bn]

                paraT.node[n]['rn'] = paraT.node[n]['rn'] - 1
                tmpp = ind0 + tmp

                if paraT.node[n]['rn'] == 0:
                    paraT.node[n]['nrft'] = 0
                    adaptFt(n, paraT, nrft)
                    for k in range(nrft):
                        try:
                            t = np.delete(t, tmpp, axis=ax)
                        except:
                            raise(ValueError)

                    Tmp = []
                    for i in range(len(forVals)):
                        ttmp = tmpp + nrft
                        if forVals[i][ax] >= tmpp:
                            if forVals[i][ax] < ttmp:
                                Tmp.append(i)
                            else:
                                forVals[i][ax] = forVals[i][ax] - nrft
                    count = 0
                    for i in Tmp:  # Tmp has to be sorted
                        forVals.pop(i-count)
                        count = count + 1

                    return(t, tmpp, nrft)

                else:

                    return(adaptDistTab(t, xx, xxt, forVals,
                                        tmpp, paraT, r, ax))

            else:

                xx = x - fbn
                xxt = lftInd[xt - bn]

                paraT.node[n]['ln'] = paraT.node[n]['ln']-1
                tmpp = ind0 + fbn

                if paraT.node[n]['ln'] == 0:

                    paraT.node[n]['nlft'] = 0
                    adaptFt(n, paraT, nlft)

                    for k in range(nlft):
                        try:
                            t = np.delete(t, tmpp, axis=ax)
                        except:
                            raise(ValueError)

                    Tmp = []
                    for i in range(len(forVals)):
                        ttmp = tmpp + nlft
                        if forVals[i][ax] >= tmpp:
                            if forVals[i][ax] < ttmp:
                                Tmp.append(i)
                            else:
                                forVals[i][ax] = forVals[i][ax] - nlft

                    count = 0
                    for i in Tmp:
                        forVals.pop(i-count)
                        count = count + 1

                    return (t, tmpp, nlft)

                else:

                    return(adaptDistTab(t, xx, xxt, forVals,
                                        tmpp, paraT, l, ax))
        # xt < tmpXT
        else:

            paraT.node[n]['fbn'] = paraT.node[n]['fbn']-1

            tt = np.delete(t, ind0 + x, axis=ax)
            adaptFt(n, paraT, 1)

            Tmp = []
            tmp = ind0 + x

            for i in range(len(forVals)):
                if forVals[i][ax] == tmp:
                    Tmp.append(i)
                elif forVals[i][ax] > tmp:
                    forVals[i][ax] = forVals[i][ax] - 1

            count = 0
            for i in Tmp:
                forVals.pop(i-count)
                count = count + 1
            # print "count", count
            # print tt,tmp,1
            return (tt, tmp, 1)

# given n (=bn) pair of transcripts to be determined
# dim[0] transcripts on one side and dim[1] transcripts on the other
# the list of forbidden cells forVal
# a graph is constructed with a source 's' and a sink 't'
# where the source is connected to dim[0] nodes
# dim[1] nodes are cnnected to the sink, and
# dim[0] nodes are connected to dim[1] nodes(all possible pairs)
# with capacity zero in case of a forbidden cell, 1 otherwise
# the function calls maximum_flow_value from networkx package
# which determines the number of paths with maximum capacity
# with the constraint that each intermediate node must appear in only one path
# returns true if the number of paths is higher or equal to n(=bn)


def maxAssign(n, dim, forVals):
    # print("max assign")
    l = len(forVals)
    if n == 0 or l < min(dim):
        return (True)
    else:
        G = nx.DiGraph()
        rows = dim[0]
        cols = dim[1]
        for i in range(rows):
            G.add_edge('s', i, capacity=1)
            for j in range(cols):
                G.add_edge(j+rows, 't', capacity=1)
                G.add_edge(i, j+rows, capacity=1)
        for f in forVals:
            G[f[0]][f[1]+rows]['capacity'] = 0
        # print "what the fuck?"
        # print G.nodes()
        # print G.edges()
        # print forVals
        # for edge in G.edges():
        #    print G[edge[0]][edge[1]]['capacity']
        try:
            # for some mysterious reason, this function
            # crashes for some particular graphs
            flowVal = nx_utils.maximum_flow_value(G, 's', 't')
        except:
            print "I could not flow"
            # print n
        # print("flow = " + str(flowVal)+" and n="+str(n))
        return(n <= flowVal)


# explore all possible solutions by branch and bound
# to find n (=bn) best pairs of transcripts
def BandB_assign(n, t, dim, pathScore, path, misRows, misCols,
                 forVals, lParaT, rParaT, lRoot, rRoot):
    global bestScore
    global bestAssign
    # print "assigning B"
    # print bestScore, bestAssign, n, pathScore, forVals
    # test if the problem can be solved
    # meaning that given the constraints(forbidden cells and dimensions)
    # we can pair n transcripts(by max flow algorithm)
    if maxAssign(n, dim, forVals):
        # if there is nothing to pair
        # then we have reduced the problem to nothing
        # we are at the end of the exploration tree
        if n == 0:
            # print "n is null"
            tmp = pathScore
            # print "tmp:",tmp,"bestScore:",bestScore
            # newPath=path[:]
            # we keep this solution if it is better than the best one so far
            if tmp < bestScore:
                bestAssign[:] = [[path, misRows, misCols]]
                bestScore = tmp
                # print "bestAssign", bestAssign
                # mR = copy.deepcopy(misRows)
                # mC = copy.deepcopy(misCols)
                # bestAssign[:] = [[newPath, mR, mC]]
            # or we add it to the list of best solutions
            elif tmp == bestScore:
                # mR = copy.deepcopy(misRows)
                # mC = copy.deepcopy(misCols)
                # bestAssign.append([newPath, mR, mC])
                bestAssign.append([path, misRows, misCols])
            # or else we do not care

        else:
            # print "n is not null"
            # find the lowest non-overlapping values in the dist array
            # x and y are the indices for the smallest value
            Mins, (x, y) = findMins(n, t, dim, forVals)
            # print Mins, (x,y)
            # the score bound is the sum of the lowest non-overlapping values
            # the best local solution along this path
            # necessarily has a higher score
            scoreBound = sum(Mins) + pathScore
            # print "scoreBound:",scoreBound, "bestScore:", bestScore
            # continue only if the score is lower than the best score
            # found so far(over all visited paths up to now)
            if scoreBound <= bestScore:
                # print "lower"
                tt = t.copy()
                newLParaT = lParaT.copy()
                newRParaT = rParaT.copy()
                newForVals = copy.deepcopy(forVals)
                newPath = path[:]
                newMisCols = copy.deepcopy(misCols)
                newMisRows = copy.deepcopy(misRows)
                newX = adaptX(x, misRows)
                newY = adaptX(y, misCols)
                # adapt the dist array for lines
                try:
                    tt, xStartDel, xNDel = adaptDistTab(
                        tt, x, newX, newForVals, 0, newLParaT, lRoot, 0)
                except:
                    print "could not update the transcript table"
                # adapt the dist array for columns
                try:
                    tt, yStartDel, yNDel = adaptDistTab(
                        tt, y, newY, newForVals, 0, newRParaT, rRoot, 1)
                except:
                    print "could not update the transcript table"
                newX = adaptX(x, misRows)  # redondant? a verifier /!\
                newY = adaptX(y, misCols)
                newPath.append((newX, newY))
                adaptMisLines(newMisRows, (xStartDel, xNDel))
                adaptMisLines(newMisCols, (yStartDel, yNDel))
                newPathScore = pathScore + min(Mins)
                newDim = (dim[0]-xNDel, dim[1] - yNDel)
                # go along this way: the problem is reduced since
                # line(s)&column(s) were removed
                # print "BandB for n-1"
                BandB_assign(n-1, tt, newDim, newPathScore, newPath,
                             newMisRows, newMisCols, newForVals,
                             newLParaT, newRParaT, lRoot, rRoot)
                # forbid the cell (x,y) and look for the second minimum
                forVals.append([x, y])
                # print "BandB for n"
                BandB_assign(n, t, dim, pathScore, path, misRows, misCols,
                             forVals, lParaT, rParaT, lRoot, rRoot)


# cette partie est ajoutee pour eviter la redondance de certains calculs
# et pour fixer un nombre de bn min pour certains noeuds en remontant les 0

# g l'arbre de gene
def leafScoreTabs(g, est, costMat, AllExons):
    f = leaves(g)
    ens = parents(g, [f, f])
    res = {}
    for e in ens[0]:
        l, r = g.successors(e)
        exSe_r = exonState(g.node[r]['trans'], AllExons)
        exSe_l = exonState(g.node[l]['trans'], AllExons)
        ll = toScore(g.node[l]['trans'], exSe_l, AllExons)
        rr = toScore(g.node[r]['trans'], exSe_r, AllExons)
        tmp = []
        exSe = est.node[e]['est']

        for i in ll:
            tmpp = []
            for j in rr:
                tmpp.append(transDist(i, j, exSe, exSe_l,
                                      exSe_r, costMat, AllExons))
            tmp.append(tmpp)
        res[str(e)] = np.array(tmp)
    return(res)

# La fonction max0 prend un tableau t et renvoie le max de 0 qui n'ont ni
# la meme ligne ni la meme colonne renvoie egalement la liste de toutes
# les combinaisons maximales de 0 de lignes et colonnes differentes

# l la liste des indices des 0 d'un tableau


def max0Rec(l, res):
    global AllPaths
    # print "res is:"
    # print res
    while len(l) > 0:
        tmp = tuple(l.pop())
        enrichPath(res, tmp)
        max0Rec(l, res)


def max0(t, res):
    l = list(np.transpose(np.where(t == 0)))
    # print l
    tmp = [[]]
    # print "tmp is:"
    # print tmp
    max0Rec(l, tmp)
    tmpRes = len(res[0])
    for i in tmp:
        len_i = len(i)
        if len_i > tmpRes:
            tmpRes = len_i
            res[:] = [i]
        elif len_i == tmpRes and len_i > 0:
            res.append(i)


# prend des listes d'indices de 0 et un nouveau couple d'indice ...
def enrichPath(l, ind):
    for i in l[:]:
        if len(i) == 0:
            i.append(ind)
        else:
            tmp = []
            count = 0
            k = 0
            for j in i:
                if j[0] != ind[0] and j[1] != ind[1]:
                    tmp.append(k)
                    count = count + 1
                k = k+1
            if count == len(i):
                i.append(ind)
            else:
                Tmp = [i[n][:] for n in tmp]
                Tmp.append(ind)
                l.append(Tmp)

# g arbre de gene, est arbre d'etats des exons,
# remonte les 0 le plus haut possible a partir des feuilles


def remonte0(g, est, costMat, AllExons):
    f = leaves(g)
    res = nx.DiGraph()
    res.add_edges_from(g.edges())
    for e in f:
        est.node[e]['est'] = exonState(g.node[e]['trans'], AllExons)
        res.node[e]['trans'] = [
            toScore(g.node[e]['trans'], est.node[e]['est'], AllExons)]
        res.node[e]['minBN'] = len(g.node[e]['trans'])
        # print e,res.node[e]['minBN']

    ens = parents(g, [f, f])
    while ens[0] != set([]):
        for e in ens[0]:
            # print e,ens[0]
            l, r = g.successors(e)
            exSe = est.node[e]['est']
            exSe_l = est.node[l]['est']
            exSe_r = est.node[r]['est']
            tmp = [[]]
            trans = [[]]
            if res.node[l]['minBN'] > 0 and res.node[r]['minBN'] > 0:
                # print "both positive"
                for i in res.node[l]['trans']:
                    for j in res.node[r]['trans']:
                        tmpp = []
                        for ii in i:
                            tmppp = []
                            for jj in j:
                                tmppp.append(transDist(ii, jj, exSe, exSe_l,
                                                       exSe_r, costMat,
                                                       AllExons))
                            tmpp.append(tmppp)

                        # max number of zeros on different lines or columns
                        # print tmpp
                        # print np.array(tmpp)
                        max0(np.array(tmpp), tmp)
                        # print tmp
                        # print trans
                        len_tmp = len(tmp[0])
                        len_trans = len(trans[0])
                        if len_tmp > len_trans:
                            for k in tmp:
                                # trans = [transFromInd(kk,i,j, exSe, exSe_l,
                                # exSe_r,costMat, AllExons) for kk in k]
                                trans = [transFromInd(k, i, j, exSe, exSe_l,
                                                      exSe_r, costMat,
                                                      AllExons)]
                        elif len_tmp == len_trans and len_tmp > 0:
                            for k in tmp:
                                TMP = transFromInd(
                                    k, i, j, exSe, exSe_l, exSe_r, costMat,
                                    AllExons)
                                trans.append(TMP)
                # print trans
                res.node[e]['trans'] = trans
                res.node[e]['minBN'] = len(trans[0])
                if e == 1000:
                    print "e values 1000"
            else:
                res.node[e]['trans'] = [[]]
                res.node[e]['minBN'] = 0
        ens = parents(g, ens)
    return(res)


def transFromInd(k, i, j, ex, ex_g, ex_d, costMat, AllExons):
    res = []
    for s in k:
        res.append(transScore(i[s[0]], j[s[1]], ex,
                              ex_g, ex_d, costMat, AllExons))
        # print res
    return(res)

############
#########
######


# Fonction principale qui affecte bn transcrits de gg a bn de dd
# bestAffectHigh pour les noeuds de hauteur au moins 2
# bestAffectLow pour les noeuds de hauteur 1
# (pour lesquels pas besoin de recalculer les couts)
def bestAffectHigh(bn, gg, dd, exSe, exSe_g, exSe_d, lpt, rpt, lRoot,
                   rRoot, costMat, AllExons):
    global bestAssign
    bestAssign = []
    global bestScore
    bestScore = 10000
    if bn == 0:
        return ([], 0)
    else:

        distMat = []
        misRows = []
        misCols = []

        for i in gg:
            tmp = []
            for j in dd:
                tmp.append(transDist(i, j, exSe, exSe_g,
                                     exSe_d, costMat, AllExons))
            distMat.append(tmp)

        t = np.array(distMat)
    # print("here ok going to BandB\n")
    # BandB_assign (bn, t, t.shape, 0, [], misRows, misCols, [], lpt, rpt,
    # lRoot, rRoot)
        try:
            BandB_assign(bn, t, t.shape, 0, [], misRows,
                         misCols, [], lpt, rpt, lRoot, rRoot)
        except:
            raise ValueError

        res = bestAssign[:]
        scoreRes = bestScore
        for a in res:
            aa = a[0]
            tmp = []
            for i in aa:
                tmp.append(t[i[0], i[1]])
            a.append(tmp)
        return (res, scoreRes)


def bestAffectLow(bn, t, lpt, rpt, lRoot, rRoot):
    global bestAssign
    bestAssign = []
    global bestScore
    bestScore = 10000
    # print "bn in bestAffectLow", bn
    if bn == 0:
        return ([])
    else:
        misRows = []
        misCols = []
    try:
        BandB_assign(bn, t, t.shape, 0, [], misRows,
                     misCols, [], lpt, rpt, lRoot, rRoot)
    except:
        print "could not assign band"
    res = bestAssign[:]
    for a in res:
        aa = a[0]
        tmp = []
        for i in aa:
            tmp.append(t[i[0], i[1]])
        a.append(tmp)
    return (res)


# res the tree, n node, c configuration in n , pc configuration in par(n)
# countBt the num of binary node in pc[0] of par(n), x the index range(bn),
# range(ln) or range(rn)
def propagateAssign(res, n, c, pc, xPar, X):

    bn = res.node[n]['bn']
    nlft = len(res.node[n]['configurations'][c]['lftInd'])
    tmpLen = len(res.node[n]['configurations'][c]['inherited'])

    elmt = 0
    found = False
    # xParTmp = 0

    if X >= bn:
        tmp = bn + nlft
        if X >= tmp:
            x = X - tmp
            s = "right"
        else:
            x = X - bn
            s = "left"
    else:
        s = "binary"
        x = X

    while not found and elmt < tmpLen:
        if res.node[n]['configurations'][c]['inherited'][elmt]['parConf'] == pc:
            res.node[n]['configurations'][c]['inherited'][elmt][s].append(
                (xPar, x))
            found = True
        else:
            elmt = elmt + 1
    if not found:
        elmt = 0
        parConfTmp = pc[:]
        tmp = {'parConf': parConfTmp, 'binary': [], 'left': [], 'right': []}
        tmp[s].append((xPar, x))
        res.node[n]['configurations'][c]['inherited'].append(tmp)

    if s != 'binary' and nx_utils.successors(res, n) != []:
        l, r = nx_utils.successors(res, n)
        PC = pc[:]
        PC.append(c)

        if s == 'left':
            xx = res.node[n]['configurations'][c]['lftInd'][x]
            lr = l
            conf = res.node[n]['configurations'][c]['lConf']
        else:
            xx = res.node[n]['configurations'][c]['rftInd'][x]
            lr = r
            conf = res.node[n]['configurations'][c]['rConf']

        propagateAssign(res, lr, conf, PC, xPar, xx)


# la fonction suivante pour ajuster un transcrit par rapport a un masque donne
# given the mask of the child and parent exons, the costs associated to
# the child exon are updated to feed the parent exon (bottom-up)
def adaptToMask(tt, parMask, mask, costMat):

    # C0 = costMat[0]
    C1_0 = costMat[1]
    C1_1 = costMat[2]
    C1_2 = costMat[3]
    t = []
    i = 0
    for ee in tt:
        e = ee[:]
        # if the child is absent
        if mask[i] == 0:
            # the parent exon is alternative
            if parMask[i] == 1:
                # the child exon is alternative and no
                e[1] = C1_1[0][1][0] + e[0]
            # the parent exon is constitutive
            elif parMask[i] == 2:
                # the child exon is absent(and no)
                e[2] = C1_2[0][0] + e[0]
                e[1] = 10000
                e[0] = 10000
        # if the child is alternative
        elif mask[i] == 1:
            # if the parent exon is absent
            if parMask[i] == 0:
                e[0] = min(e[0], e[1] + C1_0[1])
                e[1] = 10000
            # if the parent exon is constitutive
            elif parMask[i] == 2:
                e[2] = min(e[0] + C1_2[1][0], e[1] + C1_2[1][1])
                e[1] = 10000
                e[0] = 10000
        # if the child exon is constitutive
        else:
            # if the parent exon is absent
            if parMask[i] == 0:
                e[0] = C1_0[2] + e[2]
                e[2] = 10000
            # if the parent is alternative
            elif parMask[i] == 1:
                e[0] = C1_1[1][0][2] + e[2]
                e[1] = C1_1[1][1][2] + e[2]
                e[2] = 10000
        t.append(e)
        i = i + 1
    return(t)

# misL une liste de couple(i,j)
# i l'indice du premier element a retirer de T
# j le nombre d'element a retirer de T a partir de i


def iniFt(T, misL):
    count = 0
    misL.sort()
    for i in misL:
        for j in xrange(i[1]):
            T.pop(i[0] - count)
        count = count + i[1]


# t a topo ans est exons state tree
# distTabs la liste des tableaux de couts pour les noeuds de hauteur 1
# fonction principale qui remonte les assignations des feuilles a la racine
def leafAssign(t, est, distTabs, costMat, AllExons):
    res = nx.DiGraph()
    res.add_edges_from(t.edges())
    f = leaves(t)
    for e in f:
        res.node[e]['ln'] = 0
        res.node[e]['rn'] = 0
        res.node[e]['bn'] = t.node[e]['n']
        dicTmp = {}
        dicTmp['lConf'] = -1
        dicTmp['rConf'] = -1
        dicTmp['Lbn'] = []
        dicTmp['lft'] = []
        dicTmp['rft'] = []
        dicTmp['lftInd'] = []
        dicTmp['rftInd'] = []
        dicTmp['inherited'] = []
        dicTmp['Bt'] = toScore(t.node[e]['trans'], est.node[e]['est'],
                               AllExons)
        dicTmp['scoreBt'] = [0 for s in range(t.node[e]['n'])]
        res.node[e]['trans'] = t.node[e]['trans']
        res.node[e]["configurations"] = [dicTmp]

    ens = parents(t, [f, f])
    # print "ens",ens[0]
    # first level, distance matrices always the same
    for e in ens[0]:

        l, r = nx_utils.successors(t, e)
        bn = t.node[e]['bn']
        ln = t.node[e]['ln']
        rn = t.node[e]['rn']
        nl = t.node[l]['n']
        nr = t.node[r]['n']
        res.node[e]['ln'] = ln
        res.node[e]['rn'] = rn
        res.node[e]['bn'] = bn
        exSe = est.node[e]['est']
        exSe_g = est.node[l]['est']
        exSe_d = est.node[r]['est']
        res.node[e]['configurations'] = []
        countConf = 0
        distT = distTabs[str(e)]
        lParaT = generateParaT(res, l, 0, [])
        rParaT = generateParaT(res, r, 0, [])
        try:
            # print "I try"
            affectations = bestAffectLow(bn, distT, lParaT, rParaT, l, r)
        except:
            print "could not apply best low affectation"
        tmpG = res.node[l]['configurations'][0]['Bt']
        tmpD = res.node[r]['configurations'][0]['Bt']
        if bn == 0:
            dicTmp = {}
            dicTmp['lConf'] = 0
            dicTmp['rConf'] = 0
            dicTmp['Lbn'] = []
            dicTmp['Bt'] = []
            dicTmp['scoreBt'] = []
            dicTmp['inherited'] = []

            if ln > 0:
                dicTmp['lft'] = [adaptToMask(
                    m, exSe, est.node[l]['est'], costMat) for m in tmpG]
                dicTmp['lftInd'] = range(nl)
            else:
                dicTmp['lft'] = []
                dicTmp['lftInd'] = []
            if rn > 0:
                dicTmp['rft'] = [adaptToMask(
                    m, exSe, est.node[r]['est'], costMat) for m in tmpD]
                dicTmp['rftInd'] = range(nr)
            else:
                dicTmp['rft'] = []
                dicTmp['rftInd'] = []
            res.node[e]["configurations"].append(dicTmp)

        # print "affectations", affectations
        for a in affectations:
            dicTmp = {}
            dicTmp['lConf'] = 0
            dicTmp['rConf'] = 0
            dicTmp['Lbn'] = copy.deepcopy(a[0])
            dicTmp['Bt'] = []
            dicTmp['scoreBt'] = a[3]
            countBt = 0
            dicTmp['inherited'] = []

            parConf = [countConf]
            # print "parConf",parConf

            for k in a[0]:
                dicTmp['Bt'].append(transScore(tmpG[k[0]], tmpD[k[1]], exSe,
                                               exSe_g, exSe_d, costMat,
                                               AllExons))

                propagateAssign(res, l, 0, parConf, countBt, k[0])
                propagateAssign(res, r, 0, parConf, countBt, k[1])
                countBt = countBt + 1

            if ln > 0:
                lInd = range(nl)
                try:
                    iniFt(lInd, a[1])
                except:
                    print "iniFt is not ok"

                dicTmp['lft'] = [adaptToMask(tmpG[m], exSe, est.node[l]
                                             ['est'], costMat) for m in lInd]
                dicTmp['lftInd'] = lInd

            else:
                dicTmp['lft'] = []
                dicTmp['lftInd'] = []

            if rn > 0:
                rInd = range(nr)

                try:
                    iniFt(rInd, a[2])
                except:
                    print "iniFt is not ok"

                dicTmp['rft'] = [adaptToMask(tmpD[m], exSe, est.node[r]
                                             ['est'], costMat) for m in rInd]
                dicTmp['rftInd'] = rInd
            else:
                dicTmp['rft'] = []
                dicTmp['rftInd'] = []

            countConf = countConf + 1
            res.node[e]["configurations"].append(dicTmp)

    ens = parents(t, ens)
    while(ens[0] != set([])):

        for e in ens[0]:
            # print "hello"
            l, r = nx_utils.successors(t, e)
            bn = t.node[e]['bn']
            ln = t.node[e]['ln']
            rn = t.node[e]['rn']
            exSe = est.node[e]['est']
            exSe_g = est.node[l]['est']
            exSe_d = est.node[r]['est']
            nl = t.node[l]['n']
            nr = t.node[r]['n']
            res.node[e]['ln'] = ln
            res.node[e]['rn'] = rn
            res.node[e]['bn'] = bn
            res.node[e]['configurations'] = []
            countConf = 0
            bestConfScore = 10000
            lConf = res.node[l]['configurations']
            rConf = res.node[r]['configurations']
            propL = []
            propR = []

            for i in range(len(lConf)):

                tmpG = lConf[i]['Bt'] + lConf[i]['lft'] + lConf[i]['rft']
                lParaT = generateParaT(res, l, i, [])
                # lbn = res.node[l]['bn']
                # lnlt = len(lConf[i]['lftInd'])

                for j in range(len(rConf)):

                    tmpD = rConf[j]['Bt'] + rConf[j]['lft'] + rConf[j]['rft']
                    rParaT = generateParaT(res, r, j, [])
                    lpt = lParaT.copy()
                    # rbn = res.node[r]['bn']
                    # rnlt = len(rConf[j]['lftInd'])

                    try:
                        affectations, scoreAffect = bestAffectHigh(
                            bn, tmpG, tmpD, exSe, exSe_g, exSe_d, lpt, rParaT,
                            l, r, costMat, AllExons)
                    except:
                        print "could not apply best high affectation"

                    #####

                    TG = range(len(tmpG))
                    TD = range(len(tmpD))

                    if bn == 0:  # pour faire plus propre a mettre a part, ...
                        # ... au dessus de la deouble boucle i,j
                        bnNonNul = False
                        # print("cas zarb\n")
                        dicTmp = {}
                        dicTmp['lConf'] = i
                        dicTmp['rConf'] = j
                        dicTmp['Lbn'] = []
                        dicTmp['scoreBt'] = []
                        dicTmp['Bt'] = []
                        dicTmp['inherited'] = []

                        if ln > 0:
                            dicTmp['lft'] = [adaptToMask(
                                m, exSe, est.node[l]['est'], costMat) for m in
                                tmpG]
                            dicTmp['lftInd'] = TG[:]
                        else:
                            dicTmp['lft'] = []
                            dicTmp['lftInd'] = []
                        if rn > 0:
                            dicTmp['rft'] = [adaptToMask(
                                m, exSe, est.node[r]['est'], costMat) for m in
                                tmpD]
                            dicTmp['rftInd'] = TD[:]
                        else:
                            dicTmp['rft'] = []
                            dicTmp['rftInd'] = []
                        res.node[e]["configurations"].append(dicTmp)

                    elif scoreAffect < bestConfScore:
                        bnNonNul = True
                        countConf = 0
                        bestConfScore = scoreAffect
                        res.node[e]["configurations"] = []
                        propL = []
                        propR = []
                    elif scoreAffect == bestConfScore:
                        bnNonNul = True
                    else:
                        bnNonNul = False  # abus de language, ici bn est non
                        # ul mais on ne veut pas remonter ces cas

                    if bnNonNul:
                        for a in affectations:
                            dicTmp = {}
                            dicTmp['lConf'] = i
                            dicTmp['rConf'] = j
                            dicTmp['Lbn'] = a[0]  # copy.deepcopy(a[0])
                            dicTmp['Bt'] = []
                            dicTmp['scoreBt'] = a[3]
                            countBt = 0
                            dicTmp['inherited'] = []

                            # parConf = [countConf]

                            for k in a[0]:
                                dicTmp['Bt'].append(transScore(
                                    tmpG[k[0]], tmpD[k[1]], exSe, exSe_g,
                                    exSe_d, costMat, AllExons))
                                propL.append([i, [countConf], countBt, k[0]])
                                propR.append([j, [countConf], countBt, k[1]])
                                # propagateAssign(res,l,i,parConf,countBt,k[0])
                                # propagateAssign(res,r,j,parConf,countBt,k[1])
                                countBt = countBt + 1

                            if ln > 0:
                                lInd = TG[:]
                                try:
                                    iniFt(lInd, a[1])
                                except:
                                    print "lInd is not ok"

                                dicTmp['lft'] = [adaptToMask(
                                    tmpG[m], exSe, est.node[l]['est'],
                                    costMat) for m in lInd]
                                dicTmp['lftInd'] = lInd

                            else:
                                dicTmp['lft'] = []
                                dicTmp['lftInd'] = []

                            if rn > 0:
                                rInd = TD[:]
                                try:
                                    iniFt(rInd, a[2])
                                except:
                                    print "rInd is not ok"

                                dicTmp['rft'] = [adaptToMask(
                                    tmpD[m], exSe, est.node[r]['est'],
                                    costMat) for m in rInd]
                                dicTmp['rftInd'] = rInd

                            else:
                                dicTmp['rft'] = []
                                dicTmp['rftInd'] = []

                            countConf = countConf + 1
                            res.node[e]["configurations"].append(dicTmp)

            for h in propL:
                try:
                    propagateAssign(res, l, h[0], h[1], h[2], h[3])
                except:
                    print "could not propagate the assignment"

            for h in propR:
                try:
                    propagateAssign(res, r, h[0], h[1], h[2], h[3])
                except:
                    print "could not propagate the assignment"

        ens = parents(t, ens)
    return(res)


# the score will not be systematically computed for the entire forest
# the calculation stops when the current score is equal or greater
# then bestScore (which can be initialized to 10000 to avoid premature stops)
def tree_costRec(t, n, conf, parConf, cutTree, res, bestScore, costs):

    tmp = inherited(t, n, conf, parConf)
    bn = t.node[n]['bn'] - len(tmp[0])
    ln = t.node[n]['ln'] - len(tmp[1])
    rn = t.node[n]['rn'] - len(tmp[2])
    config = t.node[n]['configurations'][conf]
    res[0] = res[0] + (t.node[n]['ln'] + t.node[n]['rn']) * costs[1]
    cutTree.node[n]['ln'] = ln
    cutTree.node[n]['rn'] = rn
    if bn > 0:
        Bt = range(len(config['Bt']))
        ind = 0
        tmp[0].sort()
        for i in tmp[0]:
            try:
                Bt.pop(i-ind)
            except:
                print "problem removing the element"
            ind = ind + 1
        try:
            btScore = sum([config['scoreBt'][i] for i in Bt])
        except:
            print "problem computing the btScore"
        res[0] = res[0] + btScore + len(Bt) * costs[0]
    elif bn < 0:
        print("problem bn <0\n")
    if nx_utils.successors(t, n) != []:
        lParConf = parConf[:]
        rParConf = parConf[:]
        rParConf.append(conf)
        lParConf.append(conf)
        l, r = nx_utils.successors(t, n)
        rConf = config['rConf']
        lConf = config['lConf']
        if res[0] < bestScore:
            tree_costRec(t, l, lConf, lParConf, cutTree, res, bestScore, costs)
            tree_costRec(t, r, rConf, rParConf, cutTree, res, bestScore, costs)
        else:
            # print("une conf sautee au noeud " +str(n) +" \n")
            res[0] = 10000


def generateTreeRec(t, tt, conf, n):
    tt.node[n]['bn'] = t.node[n]['bn']
    tt.node[n]['ln'] = t.node[n]['ln']
    tt.node[n]['rn'] = t.node[n]['rn']
    tmp = len(t.node[n]['configurations'])
    tt.node[n]['configurations'] = [{}] * tmp
    tt.node[n]['configurations'][conf] = copy.deepcopy(
        t.node[n]['configurations'][conf])

    if nx_utils.successors(t, n) != []:
        lConf = t.node[n]['configurations'][conf]['lConf']
        rConf = t.node[n]['configurations'][conf]['rConf']
        l, r = nx_utils.successors(t, n)
        generateTreeRec(t, tt, lConf, l)
        generateTreeRec(t, tt, rConf, r)


def generateTree(t, conf):
    tt = nx.DiGraph()
    # tt.add_nodes_from(t.nodes())
    tt.add_edges_from(t.edges())
    generateTreeRec(t, tt, conf, 1)
    return (tt)


def tree_cost(t, baseScore, costs):
    # Res = baseScore
    confRes = []
    tRes = nx.DiGraph()
    cutTrees = []
    nConfs = len(t.node[1]['configurations'])
    print("number of configurations: " + str(nConfs) + " \n")
    # to avoid too large number of configurations
    # if nConfs > 10000:

    for i in range(nConfs):
        tt = generateTree(t, i)
        res = [0]
        cutTree = nx.DiGraph()
        cutTree.add_edges_from(tt.edges())
        tree_costRec(tt, 1, i, [], cutTree, res, baseScore, costs)  # root = 1
        if res[0] < baseScore:
            confRes = [i]
            tRes = tt.copy()
            baseScore = res[0]
            cutTrees = [cutTree.copy()]
        elif baseScore == res[0]:
            confRes.append(i)
            tRes = tt.copy()
            cutTrees.append(cutTree.copy())
        # else :
        #    print "une conf sautee\n"
    return(tRes, confRes, baseScore, cutTrees)

# next function takes a topology and calculates une borne inf
# nn node number ans p the number of enherited transcripts
# in the node nn


def evalTopoRec(t, nn, p, costs):
    succ = nx_utils.successors(t, nn)
    if succ == []:
        n = len(t.node[nn]['trans'])
        return (costs[0] * (n-p))
    else:
        l, r = nx_utils.successors(t, nn)
        bn = t.node[nn]['bn']
        ln = t.node[nn]['ln']
        rn = t.node[nn]['rn']
        n = bn + ln + rn
        tmpL = evalTopoRec(t, l, (ln+bn), costs)
        tmpR = evalTopoRec(t, r, (rn+bn), costs)
        # costs[1]: CD, costs[0]: CB
        # return(tmpL+tmpR+min(costs[1],1)*(ln+rn)+min(costs[0],1)*(n-p))
        return(tmpL + tmpR + costs[1] * (ln+rn) + costs[0] * (n-p))


def evalTopo(t, costs):
    return (evalTopoRec(t, 1, 0, costs))


def coherentTopo(t):
    res = True
    intNodes = set(t.nodes()) - leaves(t)
    for i in intNodes:
        l, r = nx_utils.successors(t, i)
        bn = t.node[i]['bn']
        ln = t.node[i]['ln']
        rn = t.node[i]['rn']
        nl = t.node[l]['n']
        nr = t.node[r]['n']
        if bn + ln > nl or bn + rn > nr or bn < 0 or ln < 0 or rn < 0:
            res = False
            break
    return(res)


# determine the types of elementary moves
# allowed at each node of the forest structure
def treeMoves(t, minBNTree):

    tt = nx.DiGraph()
    f = leaves(t)
    edges = nx_utils.get_edge_list(t)
    TMP = []
    for i in range(len(edges)):
        e = edges[i]
        if e[0] in f or e[1] in f:
            TMP.append(i)
    ind = 0
    for i in TMP:
        edges.pop(i-ind)
        ind = ind + 1

    tt.add_edges_from(edges)

    for nn in tt.nodes():
        l, r = nx_utils.successors(t, nn)
        Ln = t.node[l]['n']
        Rn = t.node[r]['n']
        bn = t.node[nn]['bn']
        minBN = minBNTree.node[nn]['minBN']
        ln = t.node[nn]['ln']
        rn = t.node[nn]['rn']
        n = t.node[nn]['n']

        if nx_utils.predecessors(t, nn) != []:
            par = nx_utils.predecessors(t, nn)[0]
            enf = nx_utils.successors(t, par)
            nLim = t.node[par]['bn']
            if nn == enf[0]:
                nLim = nLim + t.node[par]['ln']
            else:
                nLim = nLim + t.node[par]['rn']
            if nLim == 0:
                nLim = nLim + 1
        else:
            nLim = 1

        res = []
        if bn > minBN:
            res.append(1)
            res.append(5)
            res.append(7)
            if n > nLim:
                res.append(2)

        if rn > 0:
            if n > nLim:
                res.append(3)
            if ln + bn < Ln:
                res.append(12)

        if ln > 0:
            if n > nLim:
                res.append(4)
            if rn + bn < Rn:
                res.append(11)

        if ln + bn < Ln:
            res.append(6)
        if rn + bn < Rn:
            res.append(8)

        if rn + bn < Rn and ln + bn < Ln:
            res.append(9)

        if ln > 0 and rn > 0 and n > nLim:
            res.append(10)

        tt.node[nn]['nMoves'] = len(res)

        tt.node[nn]['moves'] = res

    return (tt)

# moves must not be empty


def moveTopo(baseTopo, nodes, treeMoves):

    t = baseTopo.copy()
    N = len(nodes)
    n = random.randint(0, N-1)
    nn = nodes[n]
    nMoves = treeMoves.node[nn]['nMoves']

    if nMoves == 1:
        step = treeMoves.node[nn]['moves'].pop(0)
        treeMoves.node[nn]['nMoves'] = 0
        nodes.pop(n)
    else:
        iStep = random.randint(0, nMoves-1)
        step = treeMoves.node[nn]['moves'].pop(iStep)
        treeMoves.node[nn]['nMoves'] = nMoves - 1

    if step == 1:
        t.node[nn]['bn'] = t.node[nn]['bn'] - 1
        t.node[nn]['ln'] = t.node[nn]['ln'] + 1
        t.node[nn]['rn'] = t.node[nn]['rn'] + 1
    elif step == 2:
        t.node[nn]['bn'] = t.node[nn]['bn'] - 1
    elif step == 3:
        t.node[nn]['rn'] = t.node[nn]['rn'] - 1
    elif step == 4:
        t.node[nn]['ln'] = t.node[nn]['ln'] - 1
    elif step == 5:
        t.node[nn]['bn'] = t.node[nn]['bn'] - 1
        t.node[nn]['ln'] = t.node[nn]['ln'] + 1
    elif step == 6:
        t.node[nn]['ln'] = t.node[nn]['ln'] + 1
    elif step == 7:
        t.node[nn]['bn'] = t.node[nn]['bn'] - 1
        t.node[nn]['rn'] = t.node[nn]['rn'] + 1
    elif step == 8:
        t.node[nn]['rn'] = t.node[nn]['rn'] + 1
    elif step == 9:
        t.node[nn]['bn'] = t.node[nn]['bn'] + 1
    elif step == 10:
        t.node[nn]['bn'] = t.node[nn]['bn'] + 1
        t.node[nn]['ln'] = t.node[nn]['ln'] - 1
        t.node[nn]['rn'] = t.node[nn]['rn'] - 1
    elif step == 11:
        t.node[nn]['bn'] = t.node[nn]['bn'] + 1
        t.node[nn]['ln'] = t.node[nn]['ln'] - 1
    else:
        t.node[nn]['bn'] = t.node[nn]['bn'] + 1
        t.node[nn]['rn'] = t.node[nn]['rn'] - 1

    t.node[nn]['n'] = t.node[nn]['bn'] + t.node[nn]['ln'] + t.node[nn]['rn']
    return (t, nn)


def elagMaxTopo(tt, e, n):
    succ = nx_utils.successors(tt, e)
    nn = 0
    if len(succ) > 0:
        l, r = succ
        if tt.node[e]['ln'] > 0:
            tmp = tt.node[e]['ln']
            if tmp >= n:
                tt.node[e]['ln'] = tmp - n
                tt.node[e]['n'] = tt.node[e]['n'] - n
                nn = n
            else:
                tt.node[e]['ln'] = 0
                tt.node[e]['n'] = tt.node[e]['n'] - tmp
                nn = tmp
            if nn > 0:
                elagMaxTopo(tt, l, nn)
        elif tt.node[e]['rn'] > 0:
            tmp = tt.node[e]['rn']
            if tmp >= n:
                tt.node[e]['rn'] = tmp - n
                tt.node[e]['n'] = tt.node[e]['n'] - n
                nn = n
            else:
                tt.node[e]['rn'] = 0
                tt.node[e]['n'] = tt.node[e]['n'] - tmp
                nn = tmp
            if nn > 0:
                elagMaxTopo(tt, r, nn)


def maxTopo(gt):
    tt = nx.DiGraph()
    tt.add_edges_from(gt.edges())
    tmp = leaves(gt)
    ens = [tmp, tmp]

    for i in tmp:
        tt.node[i]['trans'] = gt.node[i]['trans'][:]
        tt.node[i]['n'] = len(gt.node[i]['trans'])

    ens = parents(gt, ens)

    while ens[0] != set([]):
        for e in ens[0]:
            l, r = nx_utils.successors(gt, e)
            ml = tt.node[l]['n']
            mr = tt.node[r]['n']
            bn = min(ml, mr)
            tt.node[e]['bn'] = bn

            ln = ml-bn
            tt.node[e]['ln'] = ln

            rn = mr-bn
            tt.node[e]['rn'] = rn

            tt.node[e]['n'] = bn+ln+rn

        ens = parents(gt, ens)

    if tt.node[1]['ln'] > 0:
        elagMaxTopo(tt, 2, tt.node[1]['ln'])
        tt.node[1]['n'] = tt.node[1]['n'] - tt.node[1]['ln']
        tt.node[1]['ln'] = 0
    elif tt.node[1]['rn'] > 0:
        elagMaxTopo(tt, 3, tt.node[1]['rn'])
        tt.node[1]['n'] = tt.node[1]['n'] - tt.node[1]['rn']
        tt.node[1]['rn'] = 0
    # il reste des ln et bn en plus qu'on pourrait enlever pour
    # eventuellement baisser le score
    return (tt)


def setTopo(gt, nodes):
    tt = nx.DiGraph()
    tt.add_edges_from(gt.edges())
    tmp = leaves(gt)
    ens = [tmp, tmp]

    for i in tmp:
        tt.node[i]['trans'] = gt.node[i]['trans'][:]
        tt.node[i]['n'] = len(gt.node[i]['trans'])

    ens = set(gt.nodes()) - set(tmp)

    for e in ens:
        tt.node[e]['bn'] = nodes[e][0]
        tt.node[e]['ln'] = nodes[e][1]
        tt.node[e]['rn'] = nodes[e][2]
        tt.node[e]['n'] = nodes[e][0]+nodes[e][1]+nodes[e][2]
    return (tt)


# convert a list of exon cost tables s into a readable string
# problem: if there are multiple minima, this may not work!!
def scoreToTrans(s, AllExons):
    res = ""
    count = 0
    for i in s:
        tmp = argMin(i)
        if tmp > 0:
            res = res + AllExons[count]
        count = count + 1
    return(res)

# convert the a transcript's list of score tables
# into a readable string
# t: forest of trees
# c: configuration number
# i: node number
# j: transcript number
# k: subnode type number


def scoreToTrans2(t, c, i, j, k, AllExons):
    # print "tanscript "+str(k)
    keyInd = ['binary', 'left', 'right']
    keyTrans = ['Bt', 'lft', 'rft']
    n = len(t.nodes())
    inherTrans = inherited(t, i, c, [c]*n)
    dat = t.node[i]['configurations'][c]
    # get the transcript
    trans = dat[keyTrans[k]][j]
    convTrans = ""
    # for all exons
    for e in range(len(AllExons)):
        # print "Exon", AllExons[e],trans[e], argMinMul(trans[e])
        # get the index(ices) of the minimum(a)
        tmp = argMinMul(trans[e])
        # if only one minimum, no ambiguity
        if len(tmp) == 1:
            indMin = tmp[0]
        # if multiple minima, must determine which one to keep
        # the only possible case is alternative absent or present
        elif len(tmp) > 1 and set(tmp) == set([0, 1]):
            # if the transcript is inherited from the parents
            if j in inherTrans[k]:
                found = False
                p = 0
                # find the first assignment of the transcript
                # as it is inherited, we should find it!
                while not found:
                    datIn = dat['inherited'][p]
                    # print datIn
                    # visit all couples no matter what
                    for couple in datIn[keyInd[k]]:
                        # if the index of the transcript is found
                        # (second index of the couple)
                        if couple[1] == j:
                            found = True
                            parent = i
                            # go up along the phylogeny
                            for parConf in datIn['parConf']:
                                parent = nx_utils.predecessors(t, parent)[0]
                                # get the parent transcript
                            parTrans = t.node[parent]['configurations'][
                                parConf]['Bt'][couple[0]]
                            parIndMin = argMinMul(parTrans[e])
                            # if the parent transcript has one minimum
                            # take it as the minimum of the current
                            # transcript (child)
                            if len(parIndMin) == 1:
                                indMin = min(1, parIndMin[0])
                                # print "indMin", indMin
                            # else:
                                # print "WARNING The parent is also ambiguous!"
                    # increment the index of the config
                    p = p + 1
            # otherwise, we do not care, simply take the first one(absent)
            else:
                indMin = tmp[0]
        if indMin > 0:
            convTrans = convTrans + AllExons[e]
    return(convTrans)


# given a configuration c of a node and the list of all exons
# convert the different attributes of c
def confToTrans(c, AllExons):
    res = {}
    # list of binary transcripts as readable strings
    res['Bt'] = [scoreToTrans(i, AllExons) for i in c['Bt']]
    # list of left transcripts as readable strings
    res['Lt'] = [scoreToTrans(i, AllExons) for i in c['lft']]
    # list of right transcripts as readable strings
    res['Rt'] = [scoreToTrans(i, AllExons) for i in c['rft']]
    # the following attributes are simply copied
    res['lConf'] = c['lConf']
    res['rConf'] = c['rConf']
    res['LtInd'] = c['lftInd'][:]
    res['RtInd'] = c['rftInd'][:]
    res['BtInd'] = c['Lbn'][:]
    return(res)


def confToTrans2(t, i, c, AllExons):
    res = {}
    dat = t.node[i]['configurations'][c]
    keyTrans = ['Bt', 'lft', 'rft']
    # for all types of transcripts of node i
    for k in range(len(keyTrans)):
        # print keyTrans[k]
        res[keyTrans[k]] = []
        # for all transcripts of node i of type k
        for j in range(len(dat[keyTrans[k]])):
            res[keyTrans[k]].append(scoreToTrans2(t, c, i, j, k, AllExons))
    res['lConf'] = dat['lConf']
    res['rConf'] = dat['rConf']
    res['LtInd'] = dat['lftInd'][:]
    res['RtInd'] = dat['rftInd'][:]
    res['BtInd'] = dat['Lbn'][:]
    return(res)


# convert the list of exon cost tables into readable transcripts
# creates a new tree that comprises the same edges as the original one
# creates an attribute configurations that will look the same as
# the original one except that transcripts are displayed a strings
# and not a lists of exon cost tables
def transTree(t, AllExons):
    res = nx.DiGraph()
    res.add_edges_from(t.edges())
    for i in res.nodes():
        # print "I convert node"+str(i)
        # res.node[i]['configurations'] =
        #     [confToTrans(k, AllExons) for  k in t.node[i]['configurations'] ]
        res.node[i]['configurations'] = [confToTrans2(
            t, i, k, AllExons) for k in
            xrange(len(t.node[i]['configurations']))]
    # print "transTree finished."
    return(res)

    ###########################################################################

    # Amelioration du temps de calcul en combinant leaf assignment et
    # parcours des topologies voisines

    ###########################################################################

# returns the children of e , his brother and the brothers of all his ancestors


def siblings(res, e):
    tmp = nx_utils.successors(res, e)
    tmp.append(e)
    tmpNode = e
    ancestr = nx_utils.predecessors(res, tmpNode)
    while len(ancestr) > 0:
        par = ancestr[0]
        enf = nx_utils.successors(res, par)
        if enf[0] == tmpNode:
            tmp.append(enf[1])
        else:
            tmp.append(enf[0])
        tmpNode = par
        ancestr = nx_utils.predecessors(res, tmpNode)
    return(tmp)


# function that cleans a tree resulting from leafAssignRec
# cleaning the inherited attributes to start a new
# assignment from e to the root
def updateInherRec(res, e, n):
    nConfs = len(res.node[e]['configurations'])
    for c in range(nConfs):
        tmp = []
        for i in res.node[e]['configurations'][c]['inherited']:
            if len(i['parConf']) < n:
                tmp.append(copy.deepcopy(i))
        res.node[e]['configurations'][c]['inherited'] = tmp

    Tmp = nx_utils.successors(res, e)
    if len(Tmp) > 0:
        l, r = Tmp
        updateInherRec(res, l, (n+1))
        updateInherRec(res, r, (n+1))

# fonction analogue a leaf assign mais ici on suppose que les affectations
# ont deja ete remonte jusqu'a la racine, information disponbiles dans res
# pour un nouvelle topologie topo ne changeant qu'en un noeud e de l'ancienne
# ( l'ancienne topologie = topologie de res en input avant modifications )
# on change res qu'en redeterminant les affectations de e a la racine
# est est l'arbre des etats d'exons au niveau du gene


def leafAssignRec(t, res, est, e, distTabs, costMat, AllExons):

    l, r = nx_utils.successors(t, e)
    bn = t.node[e]['bn']
    ln = t.node[e]['ln']
    rn = t.node[e]['rn']
    nl = t.node[l]['n']
    nr = t.node[r]['n']
    res.node[e]['ln'] = ln
    res.node[e]['rn'] = rn
    res.node[e]['bn'] = bn
    exSe = est.node[e]['est']
    exSe_g = est.node[l]['est']
    exSe_d = est.node[r]['est']

    # le cas des noeud de niveau 1 (juste au dessus des feuilles)
    if nx_utils.successors(t, l) == [] and nx_utils.successors(t, r) == []:

        res.node[e]['configurations'] = []
        countConf = 0
        distT = distTabs[str(e)]
        lParaT = generateParaT(res, l, 0, [])
        rParaT = generateParaT(res, r, 0, [])

        affectations = bestAffectLow(bn, distT, lParaT, rParaT, l, r)

        tmpG = res.node[l]['configurations'][0]['Bt']
        tmpD = res.node[r]['configurations'][0]['Bt']

        if affectations == []:
            dicTmp = {}
            dicTmp['lConf'] = 0
            dicTmp['rConf'] = 0
            dicTmp['Lbn'] = []
            dicTmp['Bt'] = []
            dicTmp['inherited'] = []

            if ln > 0:
                dicTmp['lft'] = [adaptToMask(
                    m, exSe, est.node[l]['est'], costMat) for m in tmpG]
                dicTmp['lftInd'] = range(nl)
            else:
                dicTmp['lft'] = []
                dicTmp['lftInd'] = []
            if rn > 0:
                dicTmp['rft'] = [adaptToMask(
                    m, exSe, est.node[r]['est'], costMat) for m in tmpD]
                dicTmp['rftInd'] = range(nr)
            else:
                dicTmp['rft'] = []
                dicTmp['rftInd'] = []
                res.node[e]["configurations"].append(dicTmp)

        for a in affectations:
            dicTmp = {}
            dicTmp['lConf'] = 0
            dicTmp['rConf'] = 0
            dicTmp['Lbn'] = copy.deepcopy(a[0])
            dicTmp['Bt'] = []
            countBt = 0
            dicTmp['inherited'] = []

            parConf = [countConf]

            for k in a[0]:
                dicTmp['Bt'].append(transScore(tmpG[k[0]], tmpD[k[1]], exSe,
                                               exSe_g, exSe_d, costMat,
                                               AllExons))

                propagateAssign(res, l, 0, parConf, countBt, k[0])
                propagateAssign(res, r, 0, parConf, countBt, k[1])
                countBt = countBt + 1

            if ln > 0:
                lInd = range(nl)
                try:
                    iniFt(lInd, a[1])
                except:
                    print "lInd is not ok"

                dicTmp['lft'] = [adaptToMask(tmpG[m], exSe, est.node[l]
                                             ['est'], costMat) for m in lInd]
                dicTmp['lftInd'] = lInd

            else:
                dicTmp['lft'] = []
                dicTmp['lftInd'] = []

            if rn > 0:
                rInd = range(nr)

                try:
                    iniFt(rInd, a[2])
                except:
                    print "iniFt is not ok"

                    dicTmp['rft'] = [adaptToMask(
                        tmpD[m], exSe, est.node[r]['est'], costMat) for
                        m in rInd]
                    dicTmp['rftInd'] = rInd
            else:
                dicTmp['rft'] = []
                dicTmp['rftInd'] = []

            countConf = countConf + 1
            res.node[e]["configurations"].append(dicTmp)

    else:

        res.node[e]['configurations'] = []
        countConf = 0
        lConf = res.node[l]['configurations']
        rConf = res.node[r]['configurations']

        for i in range(len(lConf)):

            tmpG = lConf[i]['Bt'] + lConf[i]['lft'] + lConf[i]['rft']
            lParaT = generateParaT(res, l, i, [])
            # lbn = res.node[l]['bn']
            # lnlt = len(lConf[i]['lftInd'])

            for j in range(len(rConf)):

                tmpD = rConf[j]['Bt'] + rConf[j]['lft'] + rConf[j]['rft']
                rParaT = generateParaT(res, r, j, [])
                lpt = lParaT.copy()
                # rbn = res.node[r]['bn']
                # rnlt = len(rConf[j]['lftInd'])

                try:
                    affectations = bestAffectHigh(
                        bn, tmpG, tmpD, exSe, exSe_g, exSe_d, lpt, rParaT, l,
                        r, costMat, AllExons)
                except:
                    print "could not apply best high affectation"
                TG = range(len(tmpG))
                TD = range(len(tmpD))

                if affectations == []:
                    dicTmp = {}
                    dicTmp['lConf'] = i
                    dicTmp['rConf'] = j
                    dicTmp['Lbn'] = []
                    dicTmp['Bt'] = []
                    dicTmp['inherited'] = []

                    if ln > 0:
                        dicTmp['lft'] = [adaptToMask(
                            m, exSe, est.node[l]['est'], costMat) for
                            m in tmpG]
                        dicTmp['lftInd'] = TG[:]
                    else:
                        dicTmp['lft'] = []
                        dicTmp['lftInd'] = []
                    if rn > 0:
                        dicTmp['rft'] = [adaptToMask(
                            m, exSe, est.node[r]['est'], costMat) for
                            m in tmpD]
                        dicTmp['rftInd'] = TD[:]
                    else:
                        dicTmp['rft'] = []
                        dicTmp['rftInd'] = []
                    res.node[e]["configurations"].append(dicTmp)

                for a in affectations:
                    dicTmp = {}
                    dicTmp['lConf'] = i
                    dicTmp['rConf'] = j
                    dicTmp['Lbn'] = copy.deepcopy(a[0])
                    dicTmp['Bt'] = []
                    countBt = 0
                    dicTmp['inherited'] = []

                    parConf = [countConf]

                    for k in a[0]:
                        dicTmp['Bt'].append(transScore(tmpG[k[0]], tmpD[k[1]],
                                                       exSe, exSe_g, exSe_d,
                                                       costMat, AllExons))

                        propagateAssign(res, l, i, parConf, countBt, k[0])
                        propagateAssign(res, r, j, parConf, countBt, k[1])
                        countBt = countBt + 1

                    if ln > 0:
                        lInd = TG[:]
                        try:
                            iniFt(lInd, a[1])
                        except:
                            print "iniFt is not ok"

                        dicTmp['lft'] = [adaptToMask(
                            tmpG[m], exSe, est.node[l]['est'], costMat) for
                            m in lInd]
                        dicTmp['lftInd'] = lInd

                    else:
                        dicTmp['lft'] = []
                        dicTmp['lftInd'] = []

                    if rn > 0:
                        rInd = TD[:]

                        try:
                            iniFt(rInd, a[2])
                        except:
                            print "iniFt is not ok"

                        dicTmp['rft'] = [adaptToMask(
                            tmpD[m], exSe, est.node[r]['est'], costMat) for m in rInd]
                        dicTmp['rftInd'] = rInd
                    else:
                        dicTmp['rft'] = []
                        dicTmp['rftInd'] = []

                    countConf = countConf + 1
                    res.node[e]["configurations"].append(dicTmp)

    tmp = nx_utils.predecessors(t, e)
    if tmp != []:
        par = tmp[0]
        leafAssignRec(t, res, est, par, distTabs, costMat, AllExons)

# get an integer that represents the topology


def getCodeTopo(topo):

    code = ""
    for uu in topo.nodes():
        if nx_utils.successors(topo, uu) != []:
            code = code + str(topo.node[uu]['bn']) + \
                str(topo.node[uu]['ln'])+str(topo.node[uu]['rn'])
    return int(code)


def mkdir_subfolder(folder_path, subfolder):
    """
    Create a subfolder inside folder_path if it doesn't exists.
    """
    full_path = os.path.join(folder_path, subfolder)
    if os.path.isdir(full_path):
        warnings.warn(full_path + " already exist.", RuntimeWarning)
    else:
        os.makedirs(full_path)
    return full_path

# (best) main function with gt a genetree with transcripts at the leaves


def bestTopology(gt, AllExons, nbIt, costs, costMat, priority, SUFF, initBest,
                 slowMode, topoStart, withMemory, outputDir):
    t0 = time.time()
    nodes = []
    borInf = 0
    est = ex_state(exState(gt, costMat, AllExons), priority)
    # we start with the biggest possible forest(max number of bn)
    coupes = 0
    minBNTree = remonte0(gt, est, costMat, AllExons)
    # print "minBNTree",minBNTree.nodes()
    # for n in minBNTree.nodes() :
    #     print str(n)+": "+str(minBNTree.node[n]["minBN"])
    distTabs = leafScoreTabs(gt, est, costMat,  AllExons)
    bestTopo = nx.DiGraph()
    cutTrees = []
    countBest = 0
    countTrees = 0
    i = 0
    go = True
    f = open(os.path.join(outputDir, 'treeSearch'+SUFF+'.txt'), 'w')

    # if no minimum cost is given (default case)
    # the search starts from the topology containing the maximum
    # number of binary nodes, i.e. minimum number of trees
    if initBest == 0:
        baseTopo = maxTopo(gt)
        bestTopo = baseTopo

        i = i+1
        print(">>> " + str(i) + "\n")
        print("leaf assign for max topo\n")
        baseTree = leafAssign(baseTopo, est, distTabs, costMat, AllExons)

        # print("tree cost for max topo\n")
        tRes, bestConf, bestSc, cutTrees = tree_cost(baseTree, 10000, costs)
        baseScore = bestSc
        f.write(str(baseScore) + " " + str(bestSc) + "\n")
        print("tree cost for max topo "+str(baseScore))

        tM = treeMoves(baseTopo, minBNTree)

        for n in tM.nodes():
            if tM.node[n]['nMoves'] > 0:
                nodes.append(n)
        if nodes == []:
            print("\n~}~}~}~}~} no neighbors")
            go = False

    # otherwise the search starts from a random jump
    else:
        bestTopos_path = mkdir_subfolder(outputDir, "bestTopos")
        betterTrees_path = mkdir_subfolder(outputDir, "betterTrees")
        bestSc = initBest
        if topoStart != {}:
            baseTopo = setTopo(gt, topoStart)
            bestTopo = baseTopo

            i = i+1
            print(">>> " + str(i) + "\n")
            print("leaf assign for start topo\n")
            baseTree = leafAssign(baseTopo, est, distTabs, costMat, AllExons)

            print("tree cost for start topo\n")
            tRes, bestConf, bestSc, cutTrees = tree_cost(baseTree,
                                                         10000, costs)
            baseScore = bestSc
            print("tree cost for start topo " + str(baseScore))
            if baseScore < initBest:
                pickPrint(baseTree, os.path.join(bestTopos_path,
                                                 "topo" + str(countTrees) +
                                                 "_" + str(baseScore) + ".pk"))
                countTrees = countTrees + 1

            tM = treeMoves(baseTopo, minBNTree)

            for n in tM.nodes():
                if tM.node[n]['nMoves'] > 0:
                    nodes.append(n)
            if nodes == []:
                print("\n~}~}~}~}~} no neighbors")
                go = False

    if withMemory:
        visited = []

    # be aware that: (1) base topologies issued by random jumps
    # are not counted in the number of iterations, and(2) cut
    # neighbors of a base topology are counted
    # this implies that the number of assignments and evaluations
    # is not the same as the number of iterations
    while i < nbIt and go:
        f.flush()
        # if no base, then jump randomly
        # until a suitable base is found
        # and perform assignment for this base
        if nodes == []:
            i = i + 1
            print(">>> " + str(i) + "\n")
            print(">>>>> random jump(s) \n")
            # jumps will accumulate until a suitable base topology is found
            # they are not counted as iterations
            while borInf > bestSc or nodes == []:
                # create a new base topology
                baseTopo = aleaTopo(gt, minBNTree)
                borInf = evalTopo(baseTopo, costs)
                # print "borInf",borInf
                # generate neighbors only if the lower bound in lower
                # than best score
                if borInf <= bestSc:
                    print(">>>>> looking for neighbors... \n")
                    tM = treeMoves(baseTopo, minBNTree)

                    for n in tM.nodes():
                        if tM.node[n]['nMoves'] > 0:
                            nodes.append(n)

            # print the topology in the standard output
            # for uu in baseTopo.nodes():
                # if nx_utils.successors(baseTopo, uu) == []:
                    # uuBN = baseTopo.node[uu]['n']
                    # uuLN = 0
                    # uuRN = 0
                # else:
                    # uuBN = baseTopo.node[uu]['bn']
                    # uuLN = baseTopo.node[uu]['ln']
                    # uuRN = baseTopo.node[uu]['rn']
            # print(">>>>> node_"+str(uu)+": bn="+str(uuBN)+";
            #       ln="+str(uuLN)+"; rn="+str(uuRN)+"\n")

            # generate base forest by assignment algorithm
            baseTree = leafAssign(baseTopo, est, distTabs, costMat, AllExons)
            # print baseTree.nodes()
            # evaluate the cost of the forest
            tRes, baseConf, baseScore, tmpCutTrees = tree_cost(
                baseTree, 10000, costs)
            # print "baseScore",baseScore
            # if the cost is lower than the best solution, record it

            if baseScore <= bestSc:
                pickPrint(baseTopo, os.path.join(bestTopos_path, "topo" +
                                                 str(countBest) + "_" +
                                                 str(baseScore) + ".pk"))
                countBest = countBest + 1
            if baseScore < initBest:
                pickPrint(baseTree, os.path.join(betterTrees_path,
                                                 "tree" + str(countTrees) +
                                                 "_" + str(baseScore) + ".pk"))
                countTrees = countTrees + 1

            if baseScore < bestSc:
                print("Changing best solution\n")
                bestTopo = baseTopo
                bestSc = baseScore
                cutTrees = tmpCutTrees
            # print the cost for the base forest
            print("tree cost for base topo "+str(baseScore))
            f.write(str(baseScore) + " " + str(bestSc) + "\n")

            print(">>>>> Neighbor search ended \n")

        # be aware that this is not a else!!
        # at this point, nodes is necessarily not empty
        # evaluate a neighbor of the base topology
        t, nn = moveTopo(baseTopo, nodes, tM)
        borInf = evalTopo(t, costs)
        if withMemory:
            code = getCodeTopo(t)
            found = code in visited
        else:
            found = False
        # whether the topology is suitable or not
        # this counts as an iteration
        i = i+1
        print(">>> " + str(i) + "\n")
        # run the assignment algorithm only if the topology is suitable
        if borInf <= bestSc and not found:
            print(">>>>>> leaf assignment ...\n")
            tmpTree = leafAssign(t, est, distTabs, costMat, AllExons)
            print(">>>>>> leaf assignment ended \n")
            print(">>>>>> tree cost ...\n")
            if slowMode:
                tRes, tmpConf, tmpScore, tmpCutTrees = tree_cost(
                    tmpTree, 10000, costs)
            else:
                tRes, tmpConf, tmpScore, tmpCutTrees = tree_cost(
                    tmpTree, baseScore, costs)
            print("tree cost for current topo "+str(tmpScore))
            print(">>>>>> tree cost ended \n")
            # f.write("~~~~~~~~~ score for this topo : " +
            #     str(tmpScore) + "\n\n")

            if withMemory:
                visited.append(code)
                # print len(visited)

            # if the score of the forest is better than the base forest score,
            # change for it and generate a new neighborhood
            if tmpScore < baseScore:
                print("Changing base\n")
                baseTopo = t.copy()
                tM = treeMoves(baseTopo, minBNTree)
                nodes = []
                for n in tM.nodes():
                    if tM.node[n]['nMoves'] > 0:
                        nodes.append(n)
                baseTree = tmpTree
                baseScore = tmpScore

                # if the cost is lower than the best solution, record it
                if tmpScore <= bestSc:
                    pickPrint(t, os.path.join(bestTopos_path, "topo" +
                                              str(countBest) + "_" +
                                              str(baseScore) + ".pk"))
                    countBest = countBest + 1

                if tmpScore < bestSc:
                    print("Changing best solution\n")
                    bestTopo = t.copy()
                    bestSc = tmpScore
                    cutTrees = tmpCutTrees
            if tmpScore < initBest:
                pickPrint(tmpTree, os.path.join(betterTrees_path,
                                                "tree"+str(countTrees) + "_" +
                                                str(tmpScore) + ".pk"))
                countTrees = countTrees + 1

            # write the score to the output file
            f.write(str(tmpScore) + " " + str(bestSc) + "\n")

        # if the topology is not suitable, it is still counted as an iteration
        else:
            print("~}~}~}~}~} unsuitable topology \n")
            coupes = coupes + 1

    ti = (time.time() - t0) / 60
    # print(" execution time = " + str(ti) + " minuts\n")

    f.close()
    return(bestTopo, bestSc, coupes, ti, cutTrees)


def bestTopoRec(baseTopo, minBNTree, costs, est, distTabs, costMat, AllExons,
                bestSc, baseScore, topoRecSearch, scoreRecSearch):
    # print scoreRecSearch
    # if baseScore == 71:
        # pickPrint(baseTopo, "bestTopo71.pk")
    tM = treeMoves(baseTopo, minBNTree)

    nodes = []
    for n in tM.nodes():
        if tM.node[n]['nMoves'] > 0:
            nodes.append(n)

    bestNeigh = False

    # for all the neighbors
    while nodes != []:
        # get the neighbor
        t, nn = moveTopo(baseTopo, nodes, tM)
        borInf = evalTopo(t, costs)
        # if the topology is suitable
        if borInf <= bestSc[0]:
            # perform leaf assignment
            print(">>>>>> leaf assignment ...\n")
            tmpTree = leafAssign(t, est, distTabs, costMat, AllExons)
            print(">>>>>> leaf assignment ended \n")
            print(">>>>>> tree cost ...\n")
            tRes, tmpConf, tmpScore, tmpCutTrees = tree_cost(
                tmpTree, baseScore, costs)
            print("tree cost for current topo "+str(tmpScore))
            print(">>>>>> tree cost ended\n")
            # f.write("~~~~~~~~~ score for this topo : " +
            #     str(tmpScore) + "\n\n")

            # if the score is lower than the base score
            if tmpScore < baseScore:
                # then change for this topo as the base
                bestNeigh = True
                print("Changing base\n")
                # if the score is better than the best one, then record it
                if tmpScore < bestSc[0]:
                    bestSc[0] = tmpScore
                    # if tmpScore == 71:
                    #   pickPrint(t, "bestTopo71.pk")
                bestTopoRec(t, minBNTree, costs, est, distTabs,
                            costMat, AllExons,
                            bestSc, tmpScore, topoRecSearch, scoreRecSearch)

    # if among the nieghbors, none had a score lower than the base score
    # then to change of base occurred, then we arrived at a leaf of
    # the search tree then simply record the base topo and its score
    if not bestNeigh:
        topoRecSearch.append(baseTopo)
        scoreRecSearch.append(baseScore)

# performs a systematic exploration of the topology space from a starting
# topology the algo will visit all the neighbors of the topology
# everytime a neighbor is found with a better score than the topo,
# then the neighbor is taken as the new base,..etc recursively


def bestWideTopology(gt, AllExons, costs, costMat, priority, SUFF,
                     initBest, topoStart):
    est = ex_state(exState(gt, costMat, AllExons), priority)
    minBNTree = remonte0(gt, est, costMat, AllExons)
    distTabs = leafScoreTabs(gt, est, costMat,  AllExons)
    bestTopo = nx.DiGraph()
    cutTrees = []
    topoRecSearch = []
    scoreRecSearch = []

    # if no minimum cost is given (default case)
    # the search starts from the topology containing the maximum
    # number of binary nodes, i.e. minimum number of trees
    if initBest == 0:
        baseTopo = maxTopo(gt)

        print("leaf assign for max topo\n")
        baseTree = leafAssign(baseTopo, est, distTabs, costMat, AllExons)

        print("tree cost for max topo\n")
        tRes, bestConf, baseScore, cutTrees = tree_cost(baseTree, 10000, costs)
        print("tree cost for max topo "+str(baseScore))

        bestTopoRec(baseTopo, minBNTree, costs, est, distTabs, costMat,
                    AllExons, baseScore, baseScore, topoRecSearch,
                    scoreRecSearch)

    # otherwise the search starts from a random jump
    else:
        bestSc = [initBest]
        if topoStart != {}:
            baseTopo = setTopo(gt, topoStart)

            print("leaf assign for start topo\n")
            baseTree = leafAssign(baseTopo, est, distTabs, costMat, AllExons)

            print("tree cost for start topo\n")
            tRes, bestConf, baseScore, cutTrees = tree_cost(
                baseTree, 10000, costs)

            print("tree cost for start topo "+str(baseScore))
            bestTopoRec(baseTopo, minBNTree, costs, est, distTabs, costMat,
                        AllExons, bestSc, baseScore, topoRecSearch,
                        scoreRecSearch)

        else:
            baseTopo = aleaTopo(gt, minBNTree)
            borInf = evalTopo(baseTopo, costs)

            while borInf > bestSc[0] or nodes == []:  # CHECK : Undefined nodes
                # create a new base topology
                baseTopo = aleaTopo(gt, minBNTree)
                borInf = evalTopo(baseTopo, costs)

            print("leaf assign for start topo\n")
            baseTree = leafAssign(baseTopo, est, distTabs, costMat, AllExons)

            print("tree cost for start topo\n")
            tRes, bestConf, baseScore, cutTrees = tree_cost(
                baseTree, 10000, costs)

            print("tree cost for start topo "+str(baseScore))
            bestTopoRec(baseTopo, minBNTree, costs, est, distTabs, costMat,
                        AllExons, bestSc, baseScore, topoRecSearch,
                        scoreRecSearch)

    bestScInd = []
    bestSc = 10000
    for i in range(len(scoreRecSearch)):
        if scoreRecSearch[i] < bestSc:
            bestSc = scoreRecSearch[i]
            bestScInd = [i]
        else:
            if scoreRecSearch[i] == bestSc:
                bestScInd.append(i)

    return(topoRecSearch[bestScInd], bestSc)


def boltzmann(x, k, T):
    res = M.exp(-x/(k*T))
    return (res)


# les prochaines fonctions servent


# T un arbre, n un noeud et x le transcrit de ce neoud
# renvoie tous les transcrits aux feuilles descendants de x
def transInduits(T, n, x, c):
    if nx_utils.successors(t, n) == []:
        res = set([])
        res.add((n, x))
        return(res)
    else:
        l, r = nx_utils.successors(t, n)
        bn = T.node[n]['bn']
        conf = T.node[n]['configurations'][c]
        lnft = len(conf['lftInd'])
        lConf = conf['lConf']
        rConf = conf['rConf']
        if x < bn:
            xl = conf['Lbn'][x][0]
            xr = conf['Lbn'][x][1]
            resL = transInduits(T, l, xl, lConf)
            resR = transInduits(T, r, xr, rConf)
            return(resL | resR)

        elif x < bn + lnft:
            xl = conf['lftInd'][x-bn]
            return (transInduits(T, l, xl, lConf))

        else:
            xr = conf['lftInd'][x-bn-lnft]
            return (transInduits(T, r, xr, rConf))


def transPartRec(T, n, c, pC, res):
    Tmp = inherited(T, n, c, pC)
    bn = T.node[n]['bn']
    conf = T.node[n]['configurations'][c]
    lConf = conf['lConf']
    rConf = conf['rConf']
    BN = set(range(bn)) - set(Tmp[0])
    for i in BN:
        res.append(transInduits(T, n, i, c))

    if nx_utils.successors(t, n) != []:
        lPC = pC[:]
        lPC.append(c)
        rPC = lPC[:]

        l, r = nx_utils.successors(t, n)
        transPartRec(T, l, lConf, lPC, res)
        transPartRec(T, r, rConf, rPC, res)


def transPart(T, c):
    res = []
    transPartRec(T, 1, c, [], res)  # root = 1
    return (res)

# l liste d'indices a modifier


def compactIndList(l, X):
    n = len(l)
    if n > 0:
        Tmp = []
        maxL = max(l)
        for i in range(maxL):
            if i not in l:
                Tmp.append(i)
        Tmp.append(maxL)
        for i in range(n):
            j = 0
            while l[i] > Tmp[j]:
                j = j + 1

            l[i] = l[i] - j + X

# reucrsive function that is first called on the root of the tree
# and goes down to each leaf for a given configuration


def elagTreeRec(t, conf, n, res, ln, rn, AllExons):
    # print "this is node"+str(n)
    if nx_utils.successors(t, n) != []:
        l, r = nx_utils.successors(t, n)
        res.add_node(l)
        res.add_node(r)
        res.add_edges_from([(n, l), (n, r)])
        c = t.node[n]['configurations'][conf]
        lConf = c['lConf']
        rConf = c['rConf']
        # print lConf,rConf
        # print c
        Lbn = c['Lbn'][:]
        # print c['lftInd']
        # print c['rftInd']
        lInd = []
        rInd = []
        # numbers of bn in the children
        lbn = t.node[l]['bn']
        rbn = t.node[r]['bn']
        # trans = [scoreToTrans(i, AllExons) for i in c['Bt']]
        trans = [scoreToTrans2(t, conf, n, i, 0, AllExons)
                 for i in range(len(c['Bt']))]
        # print c['inherited']

        lLft = len(t.node[l]['configurations'][lConf]['lftInd'])
        lRft = len(t.node[l]['configurations'][lConf]['rftInd'])

        rLft = len(t.node[r]['configurations'][rConf]['lftInd'])
        rRft = len(t.node[r]['configurations'][rConf]['rftInd'])
        # fill lInd with the indices of all left subnodes of the
        # current node (number passed as parameter)
        lln, lrn = [], []
        rln, rrn = [], []

        # TMPX : les listes d'indices dans l'ordre de Lbn lftInd et rftInd de
        #     X atteints par Lbn de n
        # TMPXInd : pour chaque Lbn de n 0, 1 ou 2 sil atteint resp un Lbn un
        #     lftInd ou rftInd de X
        # C'est pour mettre a jour les Lbn pour qu'il correspondent aux
        # indices dans XftInd apres elaguage
        TmpL, TmpLInd = [[], [], []], []
        TmpR, TmpRInd = [[], [], []], []

        # in Lbn, first are listed the binary nodes, then the left nodes and
        # then the right nodes of the successors
        # the left child has index 0, the right child has index 1
        # print "Total number of subnodes: "+str(len(Lbn)+len(ln)+len(rn))
        # print "Total number of subnodes for left child: "+str(lbn+lLft+lRft)
        # print "Total number of subnodes for right child: "+str(rbn+rLft+rRft)
        for i in Lbn:
            # print i,lbn,lLft, lRft, rbn,rLft,rRft
            if i[0] >= lbn:
                # indicates a left node of the child
                if i[0] < lLft+lbn:
                    # print "left node on the left child "+str(i[0]-lbn)
                    lln.append(i[0]-lbn)
                    TmpL[1].append(i[0])
                    TmpLInd.append(1)
                    # countLln = countLln + 1
                # indicates a right node of the child
                else:
                    # print "right node on the left child "+str(i[0]-lbn-lLft)
                    lrn.append(i[0]-lbn-lLft)
                    TmpL[2].append(i[0])
                    TmpLInd.append(2)
                    # countLrn = countLrn + 1$
            # indicates a binary node of the child
            else:
                # print "binary node on the left child"
                TmpL[0].append(i[0])
                TmpLInd.append(0)

            if i[1] >= rbn:
                if i[1] < rLft+rbn:
                    # print "left node on the right child "+str(i[1]-rbn)
                    rln.append(i[1]-rbn)
                    TmpR[1].append(i[1])
                    TmpRInd.append(1)
                    # countRln = countRln + 1
                else:
                    # print "right node on the right child "+str(i[1]-rbn-rLft)
                    rrn.append(i[1]-rbn-rLft)
                    TmpR[2].append(i[1])
                    TmpRInd.append(2)
                    # countRrn = countRrn + 1
            else:
                # print "binary node on the right child"
                TmpR[0].append(i[1])
                TmpRInd.append(0)

        for i in ln:
            indexChild = c['lftInd'][i]
            if indexChild >= lbn:
                if indexChild < lLft+lbn:
                    # print "special left node on the left child " +
                    #     str(indexChild-lbn)
                    lln.append(indexChild-lbn)
                    TmpL[1].append(indexChild)
                    TmpLInd.append(1)

                else:
                    # print "special right node on the left child " +
                    #     str(indexChild-lbn-lLft)
                    lrn.append(indexChild-lbn-lLft)
                    TmpL[2].append(indexChild)
                    TmpLInd.append(2)
            else:
                # print "binary node on the left child"
                TmpL[0].append(indexChild)
                TmpLInd.append(0)

            lInd.append(indexChild)
            # trans.append(scoreToTrans(c['lft'][i], AllExons))
            trans.append(scoreToTrans2(t, conf, n, i, 1, AllExons))
            # fill rInd with the indices of all right subnodes of the
            # current node (number passed as parameter)
        for i in rn:
            indexChild = c['rftInd'][i]
            if indexChild >= rbn:
                if indexChild < rLft+rbn:
                    # print "special left node on the right child " +
                    #     str(indexChild-rbn)
                    rln.append(indexChild-rbn)
                    TmpR[1].append(indexChild)
                    TmpRInd.append(1)
                else:
                    # print "special right node on the left child " +
                    #     str(indexChild-rbn-rLft)
                    rrn.append(indexChild-rbn-rLft)
                    TmpR[2].append(indexChild)
                    TmpRInd.append(2)
            else:
                # print "binary node on the left child"
                TmpR[0].append(indexChild)
                TmpRInd.append(0)

            rInd.append(indexChild)
            # trans.append(scoreToTrans(c['rft'][i], AllExons))
            trans.append(scoreToTrans2(t, conf, n, i, 2, AllExons))
        # print ln,rn,lInd,rInd

        for i in [1, 2]:
            compactIndList(TmpL[i], lbn+(i-1)*len(TmpL[1]))

        for i in [1, 2]:
            compactIndList(TmpR[i], rbn+(i-1)*len(TmpR[1]))

        for i in xrange(len(Lbn)):
            tmp0 = TmpL[TmpLInd[i]].pop(0)
            tmp1 = TmpR[TmpRInd[i]].pop(0)
            Lbn[i] = (tmp0, tmp1)

        for i in xrange(len(lInd)):
            k = i + len(Lbn)
            tmp0 = TmpL[TmpLInd[k]].pop(0)
            lInd[i] = tmp0
        # print len(rInd)
        # print TmpRInd
        for i in xrange(len(rInd)):
            k = i+len(Lbn)
            # print i,k,TmpRInd[k]
            # print TmpR
            tmp0 = TmpR[TmpRInd[k]].pop(0)
            rInd[i] = tmp0

        # print Lbn
        # https://networkx.github.io/documentation/latest/release/migration_guide_from_1.x_to_2.0.html?
        res.node[n].update({'Lbn': Lbn, 'lftInd': lInd,
                            'rftInd': rInd, 'trans': trans})
        # print res.node[n]

        lln.sort()
        lrn.sort()
        rln.sort()
        rrn.sort()
        elagTreeRec(t, lConf, l, res, lln, lrn, AllExons)
        elagTreeRec(t, rConf, r, res, rln, rrn, AllExons)
    else:
        tmp = []
        for i in range(t.node[n]['bn']):
            tmp.append((-1, -1))
        # https://networkx.github.io/documentation/latest/release/migration_guide_from_1.x_to_2.0.html?
        res.node[n].update({'Lbn': tmp, 'lftInd': [], 'rftInd': [],
                            'trans': t.node[n].get("trans", "NA")})


def elagTree(t, conf, AllExons):
    res = nx.DiGraph()
    elagTreeRec(t, conf, 1, res, [], [], AllExons)
    return(res)
