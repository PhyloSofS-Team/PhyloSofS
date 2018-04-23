# Copyright (c) 2014-2018: Adel Ait-hamlat, Diego Javier Zea and Elodie Laine.
# This code is part of the phylosofs package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

# Modele d'inference de phylogenies de transcrits
# Pour etudier l'apparition et la fixation
# d'evenements d'epissage au cours de l'evolution

# initialisation of the data

from __future__ import division
import networkx as nx
import nx_utils

# find where to split the string
# or substring representing a binary tree or subtree
# in order to get the two children


def findSplit(treeStr):
    index = 0
    start = 0
    end = len(treeStr)
    while index == 0 and start > -1:
        i = treeStr.find(",", start, end)
        if treeStr[0:i].count("(") == treeStr[0:i].count(")") and \
            treeStr[i+1:len(treeStr)].count("(") == \
                treeStr[i+1:len(treeStr)].count(")"):
            index = i
        else:
            if i == -1:
                start = i
            else:
                start = i+1
    return(index)

# Given a string representing a tree in Newick format,
# and a list of transcripts
# Create a NetworkX graph


def convertToGraph(treeStr, transList):
    treeStr = treeStr.replace(" ", "")
    treeStr = treeStr.rstrip()
    assert treeStr[0] == "(" and treeStr[-1] == ")" and treeStr.count(
        "(") == treeStr.count(")"), "Incorrect tree"
    edgesList = []
    leafs = {}
    anc = [1]
    curr = 1
    treeStr = treeStr[1:-1]
    indexSplit = findSplit(treeStr)
    assert indexSplit > 0, "Incorrect tree"
    treeList = [treeStr[0:indexSplit], treeStr[indexSplit+1:len(treeStr)]]
    ancSuc = [treeList]
    while len(anc) > 0:
        newAnc = []
        newAncSuc = []
        for a in range(len(anc)):
            assert len(ancSuc[a]) == 2, "Incorrect tree"
            for i in range(2):
                curr = curr+1
                edgesList.append((anc[a], curr))
                if ancSuc[a][i].count(",") == 0:
                    leafs[ancSuc[a][i]] = curr
                else:
                    newAnc.append(curr)
                    tmp = ancSuc[a][i][1:-1]
                    indexSplit = findSplit(tmp)
                    assert indexSplit > 0, "Incorrect tree"
                    newAncSuc.append(
                        [tmp[0:indexSplit], tmp[indexSplit+1:len(tmp)]])
        anc = newAnc
        ancSuc = newAncSuc
    treeGraph = nx.DiGraph()
    treeGraph.add_edges_from(edgesList)
    for s in transList.keys():
        treeGraph.node[leafs[s]]['trans'] = transList[s]
    return(treeGraph)

# read input data from a text file
# that contains the geen tree in Newick format (string)
# and the list of transcripts for current species


# TO DO : Ensure correct format
def readInputDat(intree, infile):
    with open(intree, 'r') as t:
        lines = t.readlines()
        treeStr = lines[0].rstrip()

    with open(infile, 'r') as f:
        transList = {}
        for line in f:
            line = line.rstrip()
            fields = line.split(' ')
            if len(fields) >= 2:
                transList[fields[0]] = fields[1:]

    return(convertToGraph(treeStr, transList))


# get the list of all exons appearing in the leaves
def getExons(transSet):
    return [c for c in transSet.strip() if c != ':' and c != ',']

# filter based on conservation: remove exons that appear only once
# be aware that the transcripts containing these exon are not eliminated
# the number of transctips can still be reduced
# (the removed exon was the only difference)
# this is equivalent to setting a cost of zero for the any change of
# state associated to these exons


def getExonsPruned(t):
    res = {}
    transSet = getTranscripts(t)
    for tr in transSet:
        for e in tr:
            if e in res:
                res[e] += 1.0
            else:
                res[e] = 1.0
    # resPruned = [i for i in res.keys() if res.get(i) > 1]
    tbrm = [i for i in res.keys() if res.get(i) < 2]
    for n in t.nodes():
        if nx_utils.successors(t, n) == []:
            # print t.node[n]['trans']
            for k in range(len(t.node[n]['trans'])):
                for e in tbrm:
                    t.node[n]['trans'][k] = \
                        t.node[n]['trans'][k].replace(e, '')
            t.node[n]['trans'] = list(set(t.node[n]['trans']))
            # print t.node[n]['trans']

    return [i for i in res.keys() if res.get(i) > 1]

# count the number of times an exon appear


def getExonOccur(transSet):
    res = {}
    # print len(transSet)
    for t in transSet:
        for e in t:
            if e in res:
                res[e] += 1.0/len(transSet)
            else:
                res[e] = 1.0/len(transSet)
    return res

# get the leaves of the tree (indices)


def getLeaves(t):
    res = set([])
    for n in t.nodes():
        if nx_utils.successors(t, n) == []:
            res.add(n)
    return res

# get all the transcripts at the leaves


def getTranscripts(t):
    res = []
    leaves = getLeaves(t)
    for l in leaves:
        res = res + t.node[l]['trans']
    return res


def initTree(intree, infile, prune):
    Tree = readInputDat(intree, infile)
    if prune:
        exons = getExonsPruned(Tree)
    else:
        exons = getExons(getTranscripts(Tree))
    return(Tree, exons)
