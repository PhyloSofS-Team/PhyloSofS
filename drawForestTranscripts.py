# -*- coding: utf-8 -*-
"""
Created on Wed Sep  3 21:32:58 2014

A set of function to output a forest of transcripts, given a reconstructed set of trees


@author: hrichard (contributions from elaine)
"""
from sys import argv

#fichier = argv[1]

###TODO list of features:
#_ fix coloring in the exons
#_ add weblinks to ensembl
#_ get the set of

TESTING = True

import networkx as nx
import pygraphviz as pgv
import pickle
import itertools
import string
import subprocess
import os
import nx_utils

BINARY= "Lbn"
LEFT  = "lftInd"
RIGHT = "rftInd"
TOPOATTRIBUTES = [BINARY, LEFT, RIGHT]
TRANSCRIPT_DESC = ""

def ntranscripts(T, n):
    """
    Gives the number of transcripts at node n in tree T
    """
    ntrans = sum([len(T.node[n].get(x,[])) for x in TOPOATTRIBUTES])
    return ntrans

def TranscriptsTable(T,n, iconf = 0, nexons = None):
    """
    Get the transcript table of exons usage at node n of transcript T,
    for the configuration number iconf
    considering there are nexons exons at max
    """
    #List of transcripts at node (hope it is in the good order)
    lot = list(itertools.chain(*[T.node[n][x] for x in TOPOATTRIBUTES]))
    if nexons is None:
        l_ex = max(map(max, lot))
        nexons = ord(l_ex) - ord('a') + 1
    #table of transcripts usage
    tab_t = [[x if x in t else '' for x in list(string.ascii_lowercase)[:nexons]] for t in lot]
    return tab_t

def TabToRecord(tab):
    """
    output the string for a graphviz record label for the list of transcripts tab
    """
    return ["[label = \"{" + " | ".join(tr) + "}\"]" for tr in tab]


def isATree(T):
    """
    Returns true if the graph T is a tree (all nodes have at most one parent)
    """
    p = [x[1] for x in T.in_degree(T.nodes())] # in_degree & nodes instead of *_iter (NetworkX 2)
    return p.count(0) == 1 and all((x <= 1 for x in p))



def DotAllConfigurations(LT, outdir, prefix="forest_"):
    """
    From a list LT of transcripts, writes down the set of Graphviz dot files for all configurations
    in the directory outdir and compiles a pdf file with all graphs as outdir.pdf
    """
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
    for (i,T) in enumerate(LT):
        outfile="%s/%s_%d.dot" % (outdir, prefix, i)
        outpdf= "%s/%s_%d.pdf" % (outdir, prefix, i)
        ForestToDot(T, outfile)
        command = "dot -Tpdf -o %s %s" % (outpdf, outfile)
        status = subprocess.call(command, shell = True)
        ##Do something with the status
    print "Joining of all Topologies pdf is not implemented yet"


def ForestToDot(T, fileout, iconf, leafTranscripts = False, **args):
    """
    From a forest T of transcripts, writes down the Graphviz graph to fileout
    which corresponds to the configuration iconf
    Optionnally reports the transcript structure at the leaves with the parameter leafTranscripts
    Additional parameters can be passed to the graph with **args.
    """
    #TODO: get species ID for all the leaf nodes, this should be set somewhere in the tree
    #TODO: add the score of the tree where can we get this value ?
    dot_T = pgv.AGraph(directed = True, strict = True, rankdir = "TB",  newrank = True,
                       outputorder="edgesfirst", margin="0.0", splines = False, **args)
    dot_T.edge_attr.update(weight = '1', minlen = '4', dir = 'none' )
    dot_T.node_attr.update(shape = 'egg', style = 'filled', width = '1',height='0.7') # 0.1
    assert isATree(T), "Error the transcript structure provided is not a tree"
    nodes = nx.topological_sort(T)
    #we keep a dictionnary of all nodes in the forest
    forest_nodes = {}
    nspecies = []
    nbDeaths = 0
    ##Get the conf number at each node for each configuration
    ##get the total number of transcripts
    ##TODO if transcript structure is asked, add a node for each structure in the subgraph, and
    ##impose the sink level for the structure nodes under the transcript node
    fout = open(fileout+"_config"+str(iconf)+".info", "w")
    fout2 = open(fileout+"_config"+str(iconf)+".sum", "w")
    #print nodes
    for n in sorted(nodes):
        ##create as many nodes as transcripts for the given configuration
	nlist = [ "%d_%d" % (n, i) for i in range(ntranscripts(T,n)) ]
	#print T.node[n]
	for k in range(len(nlist)):
	    #print T.node[n]['trans'][k]
	    fout.write(nlist[k]+": "+T.node[n]['trans'][k]+"\n")
            dot_T.add_node(nlist[k],label="")
        fout.write("\n")
	if n>1:
	    parent = nx_utils.predecessors(T, n)[0]
	    if sorted(nx_utils.successors(T, parent)).index(n) == 0 :
	        choice = RIGHT
                #print "I'm right", n, "my parent is", parent, T.node[parent][RIGHT]
	    else:
		choice = LEFT
		#print "I'm left", n, "my parent is", parent, T.node[parent][LEFT]
            for l in T.node[parent][choice]:
	        #sl_t = "%d_%d" % (sl,l)
	        ##ajoute un noeud square pour les morts
	        ##Pbme, il est pas du bon côté
	        t_death = "%d_%d_death" % (n, l)
		nbDeaths = nbDeaths + 1
		#print t_death
	        dot_T.add_node(t_death, shape = "triangle", label = "", width='0.3')
	        nlist.append(t_death)
	        #dot_T.add_edge(dad_t, t_death)
                #i = i + 1
	        #dad_t = "%d_%d" % (n,i)
        #dot_T.add_nodes_from(nlist)
        currank = "same"
        #all species are pushed at the bottom of the tree
        if (nx_utils.get_out_degree(T, n) == 0):
            currank = "sink"
            nspecies.extend(nlist)
        dot_T.add_subgraph(nlist, name = "cluster_%d" % (n),
                           label = "", rank = currank) #"%d" % (n), rank = currank)
        for fn in nlist:
            forest_nodes[fn] = currank
    fout.close()
    ##All actual species have to be at the same level
    #dot_T.add_subgraph(nspecies, rank = "sink")
    #for fn in nspecies: forest_nodes[fn] = "sink"


    def _testcreate(s,e, d):
        ##simple logging of errors, should raise
        a = [i for (i,x) in enumerate([not d.has_key(s),not  d.has_key(e)]) if x ]
        if len(a)>0:
            print "**** Error when creating edge from %s to %s *****" % (s,e)
            for x in a: print "%s does not exist" % ([s,e][x])

    ##Drawing all the edges from the internal nodes
    for n in filter(lambda x: nx_utils.get_out_degree(T,x) > 0, nodes):
        i = 0
        for i_tr in TOPOATTRIBUTES:
            dad_t = "%d_%d" % (n, i)
            ##left and right sons (sorted left-to-right by construction)
            sl, sr = sorted(nx_utils.successors(T,n))
            if i_tr == BINARY:
                #2 sons
                for (l,r) in T.node[n][i_tr]:
                    sl_t, sr_t = "%d_%d" % (sl,l),  "%d_%d" % (sr,r)
                    dot_T.add_edges_from([(dad_t,sl_t), (dad_t, sr_t)])
                    _testcreate(dad_t,sl_t, forest_nodes)
                    _testcreate(dad_t, sr_t, forest_nodes)
                    i = i + 1
                    dad_t = "%d_%d" % (n, i)
            elif i_tr == LEFT:
                ##left sons
                for l in T.node[n][i_tr]:
                    sl_t = "%d_%d" % (sl,l)
                    dot_T.add_edges_from([(dad_t, sl_t) ])
                    _testcreate(dad_t,sl_t, forest_nodes)
                    ##ajoute un noeud square pour les morts
                    ##Pbme, il est pas du bon côté
                    t_death = "%d_%d_death" % (sr, l)
		    #print t_death
                    #dot_T.add_node(t_death, shape = "triangle", label = "", width='0.3')
                    dot_T.add_edge(dad_t, t_death)
                    i = i + 1
                    dad_t = "%d_%d" % (n,i)
            elif i_tr == RIGHT:
                ##right
                for r in T.node[n][i_tr]:
                    sr_t = "%d_%d" % (sr,r)
                    dot_T.add_edges_from([(dad_t, sr_t)])
                    _testcreate(dad_t, sr_t, forest_nodes)
                    ##ajoute un noeud square pour les morts
                    ##Pbme, il est pas du bon côté
                    t_death = "%d_%d_death" % (sl, r)
		    #print t_death
                    #dot_T.add_node(t_death, shape = "triangle", label = "", width='0.3')
                    dot_T.add_edge(dad_t, t_death)
                    i = i + 1
                    dad_t = "%d_%d" % (n,i)
            else :
                ##error
                raise "Error in drawing"
    ### TODO Still needs a function drawing all the ancestral transcripts at the
    ### nodes
    ### getting all sets of transcripts which are related (same connected component)
    ### and print them with the same color symbol
    #print "Nodes of dot_T"
    #print dot_T.nodes()
    #print "edges"
    #print dot_T.edges()

    # NetworkX 1.11: nx_agraph is not longer imported from
    ntr_totAll = len(sorted(nx.connected_components(nx.nx_agraph.from_agraph(dot_T).to_undirected())))
    ntrans = filter( lambda x: len(x) > 1,  nx.connected_components(nx.nx_agraph.from_agraph(dot_T).to_undirected()))
    ntr_tot = len(ntrans)
    fout2.write("Total number of trees: "+str(ntr_tot)+"\n")
    fout2.write("Total number of deaths: "+str(nbDeaths)+"\n")
    fout2.write("Total number of orphans: "+str(ntr_totAll-ntr_tot)+"\n")
    colorscheme = "set1%d" % (max([3, ntr_tot]))
    if ntr_tot > 9:
        colorscheme = "set3%d" % (max([3, ntr_tot]))
    Graph_args = {"colorscheme" : colorscheme} ##Additional args for graph drawing
    #orphan transcripts stay in black
    for (i,l_n) in enumerate(ntrans):
        for n_id in l_n:
            dot_T.get_node(n_id).attr.update({'colorscheme': colorscheme,'fillcolor': i+1})
    ##Options of the graph
    ###update the arguments in the graphs
    dot_T.graph_attr.update(Graph_args)
    fout = open(fileout+"_config"+str(iconf)+".dot", "w")
    source = dot_T.to_string()
    source = source.replace("\\\n", '')
    source = source.replace('''"<''', '''<''')
    source = source.replace('''>"''', '''>''')
    source = source.replace("\\\"", "\"")
    fout.write(source)
    fout.close()
    fout2.close()


#TESTING=False
#if TESTING:
#    ##Des tests temporaires avec les structures d'arbre de MAPK8
#    with open("MAPK8.transcriptTree.pickle", 'rb') as f:
#        transT = pickle.load(f)
#    with open("test_affichage.txt") as f:
#        transTmultiple = pickle.load(f)
#    with open("MAPK8.result.pickle", 'rb') as f:
#        resT = pickle.load(f)
#    with open("MAPK8.Gene.pickle", 'rb') as f:
#        MAPK8_Gene = pickle.load(f)
#
#    with open("rres.pk", 'rb') as f:
#        transT2 = pickle.load(f)
#
#    ForestToDot(transT, "test.dot")
#    status = subprocess.call("dot" + " -Tpdf -o test.pdf test.dot", shell = True)
#
#    ForestToDot(transT2[0], "tmult1.dot", 0)
#    status = subprocess.call("dot" + " -Tpdf -o tmult1.pdf tmult1.dot", shell = True)
#
#    #ForestToDot(transTmultiple, "tmult2.dot", 1)
#    #status = subprocess.call("dot" + " -Tpdf -o tmult2.pdf tmult2.dot", shell = True)

#with open(fichier, 'rb') as f:
#    transU = pickle.load(f)
#count = 0
#for i in transU :
#	ForestToDot(i, "test.dot")
#	status = subprocess.call("dot" + " -Tpdf -o Res_"+str(count)+".pdf test.dot", shell = True)
#	count = count + 1

#ForestToDot(transU[1], "test1.dot")
#status = subprocess.call("dot" + " -Tpdf -o test1.pdf test.dot", shell = True)
