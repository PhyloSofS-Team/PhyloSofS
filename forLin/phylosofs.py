#####  Modele d'inference de phylogenies de transcrits
#####  Pour etudier l'apparition et la fixation 
#####  d'evenements d'epissage au cours de l'evolution

## Par Elodie Laine, Hugues Richard, Adel Ait-hamlat
## Usage: python transPhyl3.py CB CD cm SUFF
## CB: birth cost
## CD: death cost
## cm: mutation cost
## SUFF: suffix for the output files

from __future__ import division
import sys 
import time
import random
import networkx as nx
import itertools as it
import numpy as np
import math as M
import copy
import initData	
import pickle as pk	
from inferPhylo import *
from modelIsoforms import * 


def printUsage():
    print """Usage: python phylosofs.py
             Mandatory:  
             ==========
                        -c            text file containing values for parameters (for molecular modeling part)
                        -inSeq        text file containing input data: string representing the gene tree in Newick format on the first line, 
                                      and then the list of transcripts for each leave (leaf_name: t1,t2,t3...)
                        -mode         P and/or M for phylogenetic inference and/or molecular modeling 
             Optional:
             =========
                        -b            birth cost (by default: 5)
                        -d            death cost (by default: 3)
                        -inStruct     text file containing input data: either a directory where multifasta files 
                                      are located (one file per species) or a single fasta file with the sequence
                                      of only one transcript (in that case, the uniq option must be set to TRUE)
                                      (by default, the current directory is used)
                        -m            mutation cost (by default: 2)
                        -ni           number of iterations (by default: 1)
                        -nt           number of templates to retain (by default: 5)
                        -noPrune      disable the removal of exons that appear in only one transcript
                        -o            output directory (by default: $PWD)
                        -only3D       perform only the 3D modeling step (skip template search)
                        -onlyQuality  perform only the 3D models quality assessment 
                        -printOnly    perform only the generation of the PDF file enabling to visualize a transcripts'  
                                      phylogeny given as input via the option -topo 
                        -uniq         treat only one transcript whose sequence is taken from the fasta file
                                      indicated by the -i option
                        -s            starting score (by default: not considered), if no score is given, the algorithm starts the search 
                                      by the topology corresponding to the maximum number of binary subnodes at each nodes (forest with 
                                      the smallest possible number of trees), otherwise it starts from a randomly chosen topology
                        -suff         suffix (by default: _bdm_n)
                        -topo         initial topology (by default: maximum or random topology), or transcripts'
                                      phylogeny to be printed out (if the -printOnly option is active)"""

if (__name__ == '__main__'):

    nb_arguments = len(sys.argv)
    if(nb_arguments < 2):
        printUsage()
        exit(1)

    random.seed()

    doPhylo = False
    doModel = False

    try:
        mode = sys.argv[sys.argv.index("-mode")+1]
        if "P" in mode:
            doPhylo = True
            try:
                inputFile = sys.argv[sys.argv.index("-inSeq")+1]
                os.path.isfile(inputFile)
            except:
                sys.stderr.write("You must give an existing input text file for phylogenetic inference. See usage instructions.\n\n")
                exit(2)
        if "M" in mode:
            doModel = True
            try:
                configFile = sys.argv[sys.argv.index("-c")+1]
                os.path.isfile(configFile)
                HHBLITS, ADDSS, HHMAKE, HHSEARCH, HHMODEL, PROCHECK, NACCESS, HHDB, STRUCTDB, ALLPDB, NCPU, CONTEXTLIB = init(configFile)
            except:
                sys.stderr.write("You must give an existing configuration file for molecular modeling. See usage instructions\n\n")
                exit(2)
    except:
        sys.stderr.write("You must give at least one action to be performed by PhyloSofS. See usage instructions.\n\n")
        exit(2)

    # birth cost
    try :
        CB = int(sys.argv[sys.argv.index("-b") + 1])
    except :
        CB = 5

    # death cost
    try :
        CD = int(sys.argv[sys.argv.index("-d") + 1])
    except :
        CD = 3

    try:
        pathTransSeqs = sys.argv[sys.argv.index("-inStruct")+1]
    except:
        pathTransSeqs = "./"

    # mutation cost
    try :
        cm = int(sys.argv[sys.argv.index("-m") + 1])
    except :
        cm = 2

    # number of iterations
    try :
        nbIt = int(sys.argv[sys.argv.index("-ni") + 1])
    except :
        nbIt = 1

    try:
        nbTemp = int(sys.argv[sys.argv.index("-nt")+1])
    except:
        nbTemp = 5

    # output directory
    try:
        outputDir = os.path.expanduser(sys.argv[sys.argv.index("-o")+1])
    except:
        outputDir = "./"

    # only 3D model reconstruction for modeling part
    try:
        sys.argv[sys.argv.index("-only3D")]
        only3D = True
    except :
        only3D = False

    # only assess the quality of the 3D models for modeling part
    try:
        sys.argv[sys.argv.index("-onlyQuality")]
        onlyQuality = True
    except :
        onlyQuality = False

    # initial score
    try:
        initBest = int(sys.argv[sys.argv.index("-s")+1])
        slowMode = True 
    except:
        initBest = 0
        slowMode = False

    # suffix
    try:
        SUFF = sys.argv[sys.argv.index("-suff")+1]
    except:
        SUFF = "_"+str(CB)+str(CD)+str(cm)+"_"+str(nbIt)

    # input topology
    try:
        f = sys.argv[sys.argv.index("-topo")+1]
        topoStart = pk.load(open(f,'rb'))
    except:
        topoStart = {}

    # only one transcript to be modeled
    try:
        sys.argv[sys.argv.index("-uniq")]
        uniq = True
        if '.fa' not in pathTransSeqs:
            sys.stderr.write('You must give a fasta input file with option -uniq.')
            exit(2)
    except:
        uniq = False

    # print only the solution corresponding to the input topology
    try :
        sys.argv[sys.argv.index("-printOnly")]
        printOnly = True
    except :
        printOnly = False

    try :
        sys.argv[sys.argv.index("-withMemory")]
        withMemory = True
    except :
        withMemory = False

    try :
	sys.argv[sys.argv.index("-noPrune")]
	prune = False
    except :
	prune = True

    ###################### phylogenetic inference ######################

    if doPhylo:
        
        print "--------------------------------------------"
        print "Running phylogenetic reconstruction step..."
        print "--------------------------------------------"

        bestAssign = []
        bestScore = 10000

        # couts de changement d'etat des exons pour les especes ancestrales:
        # costs are ok
        # C1 : niveau proteique (etape 2)
        #C1= [ [0,cm,0], [cm,0,0], [0,cm,0] ]
        C1_1 = [ [[0,cm,0],[cm,0,0]] , [[0,0,cm],[0,0,0]] ]
          ##   1 -> 1                 0/2 -> 1     
        C1_2 = [ [0,0,0], [cm,0,0] ]
          ##   0->2     1/2->2
        C1_0 = [0,cm,0]

        # C0 : niveau genique (etape 1)
        C0 = [ [0,1,1], [1,0,1], [1,1,0] ]

        priority = (1,0,2)

        costs = (CB, CD, cm)
        costMat=(C0,C1_0,C1_1,C1_2)

        dat,AllExons=initData.initTree(inputFile, prune)
        print "The exons are:"
        print AllExons
        nExons = range(len(AllExons))
        print nExons

        if topoStart == {}:
            res = bestTopology(dat, AllExons, nbIt, costs, costMat, priority, SUFF, initBest, slowMode, topoStart, withMemory, outputDir)
        else:
            if printOnly:
                res = [topoStart]
            else:
                res = bestWideTopology(dat, AllExons, costs, costMat, priority, SUFF, initBest,topoStart)

        if len(res[0].nodes())>0:
            exSt = ex_state(exState(dat, costMat, AllExons), priority)
            distTabs = leafScoreTabs(dat, exSt, costMat, AllExons)
            if printOnly:
                tree = topoStart
                writeOutput((tree,0,0), exSt, SUFF, costs, AllExons, outputDir)
            else:
                tree = leafAssign(res[0], exSt, distTabs, costMat, AllExons)
                writeOutput((tree,res[2],res[3]), exSt, SUFF, costs, AllExons, outputDir)
        else:
            print "No suitable topology could be found."

    ###################### molecular modeling ######################
    if doModel:

        print "--------------------------------------------"
        print "Running molecular modeling step..."
        print "--------------------------------------------"

        os.chdir(outputDir)
        # create as many directories as fasta input files
        # and inside as many fasta files as transcripts
        if not uniq and not onlyQuality:
            print 'prepare intputs...'
            prepareInputs(pathTransSeqs)

        # determine the number of templates that will be retained
        selTemp = ""
        for i in range(nbTemp):
            selTemp += str(i+1)
            selTemp += " "

        # perform the whole 3D modelling process
        dirs = [x[0] for x in os.walk('.')]
        if not onlyQuality:
            print 'launch the 3D modelling process'
            if uniq:
                mydir, trans = os.path.split(pathTransSeqs)
                os.chdir(mydir)
                runModelProcess(HHBLITS,ADDSS,HHMAKE,HHSEARCH,HHMODEL,HHDB,STRUCTDB,ALLPDB,NCPU,'./'+trans,selTemp,only3D)
                os.chdir('..')
            else:
                print dirs
                for mydir in dirs[1:] : 
                    os.chdir(mydir)
                    for trans in gl.glob('./*.fa'): 
                        runModelProcess(HHBLITS,ADDSS,HHMAKE,HHSEARCH,HHMODEL,HHDB,STRUCTDB,ALLPDB,NCPU,trans,selTemp,only3D)
                    os.chdir('..')

        # assess the quality of the models
        print 'assess the quality of the models'
        fOUT = open('quality.sum','w')
        fOUT.write('# procheck: Ideally, scores should be above -0.5. Values below -1.0 may need investigation.\n')
        fOUT.write('# dope: This is a Z-score; positive scores are likely to be poor models, while scores lower than -1 or so are likely to be native-like.\n')
        fOUT.write('transcript\tlenFull\tpercentSS\tlenModel\tdihedrals\tcovalent\toverall\tdope\trSurf\trHydroph\n')
        for mydir in dirs[1:]: 
            try:
                os.chdir(mydir)
                try:
                    os.stat('procheck/')
                except:
                    os.mkdir('procheck/')
                os.chdir('procheck/')
                for prot in gl.glob('../*.B99990001.pdb'):
                    mydir, trans = os.path.split(prot)
                    pref = trans.split('.')[0]
                    lenModel = computeLenModel('../'+pref)
                    lenFull, percentSS = computeSS('../'+pref)
                    dihedrals, covalent, overall = assessQuality(PROCHECK,prot,pref)
                    zscore = assessNormalizedDopeScore(prot)
                    rSurf, rHydroph = computeRatioSASA(NACCESS,prot,pref)
                    fOUT.write(pref + '\t' + str(lenFull) + '\t' + str(percentSS) + '\t' + str(lenModel) + '\t' + dihedrals + '\t' + covalent + '\t' + overall + '\t' + str(zscore) + '\t' + str(rSurf) + '\t' + str(rHydroph) + '\n')
                os.chdir('..')
                os.system('rm -r procheck/')
                os.chdir('..')
            except:
                print 'Problem with ' + mydir
        fOUT.close()
