# -*- coding: utf-8 -*-

# Copyright (c) 2014-2018: Elodie Laine, Hugues Richard, Adel Ait-hamlat
# and Diego Javier Zea.
# This code is part of the phylosofs package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

# Modele d'inference de phylogenies de transcripts
# Pour etudier l'apparition et la fixation
# d'evenements d'epissage au cours de l'evolution

from __future__ import division

import sys
import os
import glob
import random
import time
import pickle as pk
import argparse
import pathlib
import initData
import utils
import inferPhylo as ip
import modelIsoforms as mi
# CB: birth cost
# CD: death cost
# cm: mutation cost
# SUFF: suffix for the output files


def check_argument_groups(parser, arg_dict, group, argument, required):
    """
    Check for use of arguments.

    Raise an error if the parser uses an argument that belongs to an argument
    group (i.e. mode) but the group flag is not used or if the argument is
    required in the argument group and it's not used.
    Notes:
    - argument and group should be strings and should start with - or --
    - group should be a flag with action 'store_true'
    - argument shouldn't be a flag, and should be None if it isn't used.

    >>> import argparse
    >>> parser = argparse.ArgumentParser()
    >>> parser.add_argument('--phylo', action='store_true')
    _StoreTrueAction(...)
    >>> parser.add_argument('--inseq')
    _StoreAction(...)
    >>> args = parser.parse_args(['--phylo', '--inseq', 'not_none'])
    >>> check_argument_groups(parser, vars(args), '--phylo', '--inseq', True)
    """
    arg_name = argument.replace("-", "")
    group_name = group.replace("-", "")
    if not arg_dict[group_name]:
        if arg_dict[arg_name] is not None:
            parser.error("phylosofs requires " + group + " if " + argument + " is used.")
    else:
        if required and arg_dict[arg_name] is None:
            parser.error("phylosofs requires " + argument + " if " + group + " is used.")
    return None


def parse_command_line():
    """
    Parse command line.

    It uses argparse to parse phylosofs' command line arguments and returns the
    argparse parser.
    """
    parser = argparse.ArgumentParser(
        prog="phylosofs",
        description="""
        PhyloSofS (PHYLOgenies of Splicing isOforms Structures) is a tool
        to model the evolution and structural impact of alternative
        splicing events.
        """,
        epilog="""
        It has been developed at LCQB (Laboratory of Computational and
        Quantitative Biology), UMR 7238 CNRS, Sorbonne Université.
        If you use it, please cite:
        Ait-hamlat A, Polit L, Richard H, Laine E. Transcripts evolutionary
        conservation and structural dynamics give insights into the role of
        alternative splicing for the JNK family. bioRxiv. 2017 Jan 1:119891.
        """,
    )

    phylo_args = parser.add_argument_group(
        'Phylogenetic reconstruction', """
        Arguments used for reconstructing transcripts' phylogenies.
        """)

    model_args = parser.add_argument_group(
        'Molecular modeling', """
        Arguments used for modeling the tertiary structures of the isoforms.
        """)

    phylo_args.add_argument('-P', '--phylo',
                            help='do the phylogenetic inference',
                            action='store_true')
    model_args.add_argument('-M', '--model',
                            help='do the molecular modelling',
                            action='store_true')
    phylo_args.add_argument('-n', '--tree',
                            help='string representing the gene tree in Newick format')
    phylo_args.add_argument('-t', '--transcripts',
                            help='text file containing the list of transcripts for '
                            'each leave (leaf_name: t1,t2,t3...)')
    phylo_args.add_argument('-b',help='birth cost',type=int,default=5)
    phylo_args.add_argument('-d',help='death cost',type=int,default=3)
    model_args.add_argument('-i', '--instruct',
                            help='text file containing input data: either a directory where '
                            'multifasta files are located (one file per species) or a single '
                            'fasta file with the sequence of only one transcript (in that '
                            'case, the unique option must be set to TRUE)',
                            default='')
    phylo_args.add_argument('-m',help='mutation cost',type=int,default=2)
    phylo_args.add_argument('--ni',help='number of iterations',type=int,default=1)
    model_args.add_argument('--nt',help='number of templates to retain',type=int,default=10)
    phylo_args.add_argument('--noprune',
                            help='disable the removal of exons that appear in only one transcript',
                            action='store_true')
    parser.add_argument('-o', '--outputdir',help='output directory',default='.')
    model_args.add_argument('--only3D',
                            help='perform only the 3D modeling step (skip template search)',
                            action='store_true')
    model_args.add_argument('--onlyquality',
                            help='perform only the 3D models quality assessment',
                            action='store_true')
    phylo_args.add_argument('--printonly',
                            help='perform only the generation of the PDF file enabling to '
                            'visualize the transcripts',
                            action='store_true')
    model_args.add_argument('--uniq',
                            help='treat only one transcript whose sequence is taken from '
                            'the fasta file indicated by the -i option',
                            action='store_true')
    phylo_args.add_argument('-s',
                            help='starting score (by default: not considered), '
                            'if no score is given, the algorithm starts the search '
                            'by the topology corresponding to the maximum number of '
                            'binary subnodes at each nodes (forest with the smallest '
                            'possible number of trees), otherwise it starts from a '
                            'randomly chosen topology')
    phylo_args.add_argument('--suffix',
                            help='suffix, it is _(birth)(death)(mutation)_(iterations) by '
                            'default, e.g. _532_1')
    phylo_args.add_argument('--topo',
                            help="initial topology (by default: maximum or random topology), "
                            "or transcripts' phylogeny to be printed out "
                            "(if the -printonly option is active)")
    phylo_args.add_argument('--withmemory',
                            action='store_true')
    model_args.add_argument(
                            '--hhlib',
                            help='Path to the hh-suite',
                            default='')
    model_args.add_argument(
                            '--hhdb',
                            help="Path to the sequence database for the HH-suite, "
                            "e.g. Uniclust30",
                            default='')
    model_args.add_argument(
                            '--structdb',
                            help='Path to the structure database for the HH-suite, e.g. PDB',
                            default='')
    model_args.add_argument('--ncpu',
                            help='Number of CPUs',
                            default='1')
    model_args.add_argument('--allpdb',
                            help='Path to all the pdb database (in cif format)',
                            default='allpdb')

    args = parser.parse_args()

    if args.phylo is None and args.model is None:
        parser.error("phylosofs requires --phylo or/and --model flags.")

    arg_dict = vars(args)

    check_argument_groups(parser, arg_dict, '--phylo', '--transcripts', True)
    check_argument_groups(parser, arg_dict, '--phylo', '--tree', True)
    #check_argument_groups(parser, arg_dict, '--model', '-c', True)
    # Check flag arguments
    if not args.model:
        if args.only3D:
            parser.error("phylosofs requires --model if --only3D is used.")

    if args.model:
        if not args.instruct:
            parser.error("phylosofs requires an input fasta file if --model is used.")

    return args



def choose_path(hhlib, program):
    """
    Return the first executable path between hhlib/build/bin/program and
    hhlib/bin/program
    """
    paths = [os.path.join(hhlib, "build", "bin", program),
             os.path.join(hhlib, "bin", program)]
    for path in paths:
        if os.path.isfile(path) and os.access(path, os.X_OK):
            return path
    raise Exception("hhlib is %s, but %s aren't executable paths for "
                    "program %s" % (hhlib, paths, program))


def doit(doPhylo,
         doModel,
         inputTree,
         inputFile,
         CB,  # birth cost
         CD,  # death cost
         pathTransSeqs,
         cm,  # mutation cost
         nbIt,  # number of iterations
         nbTemp,
         outputDir,
         only3D,
         onlyQuality,
         starting_score,
         SUFF,
         topology_file,
         uniq,
         printOnly,
         withMemory,
         prune,
         HHLIB,
         HHDB,
         STRUCTDB,
         ALLPDB,
         NCPU,
         ):
    """
    doit(args.phylo,
         args.model,
         args.tree,
         args.transcripts,
         args.b,
         args.d,
         args.instruct,
         args.m,
         args.ni,
         args.nt,
         args.outputdir,
         args.only3D,
         args.onlyquality,
         args.s,
         args.suffix,
         args.topo,
         args.uniq,
         args.printonly,
         args.withmemory,
         not args.noprune,
         args.hhlib,
         args.procheck,
         args.naccess,
         args.hhdb,
         args.structdb,
         args.allpdb,
         args.ncpu,
         )
    """
    random.seed()

    outputDir = os.path.abspath(outputDir)
    pathTransSeqs = os.path.abspath(pathTransSeqs)

    if doPhylo:
        for f in [inputFile, inputTree]:
            if not os.path.isfile(f):
                sys.stderr.write(f + " : You must give a valid input "
                                 "file for phylogenetic inference.")
    if starting_score is not None:
        initBest = starting_score
        slowMode = True
    else:
        initBest = 0
        slowMode = False

    if SUFF is None:
        SUFF = "_" + str(CB) + str(CD) + str(cm) + "_" + str(nbIt)

    # input topology
    if topology_file is not None:
        topoStart = pk.load(open(topology_file, 'rb'))
    else:
        topoStart = {}

    # only one transcript to be modeled
    if uniq and '.fa' not in pathTransSeqs:
        sys.stderr.write('You must give a fasta input file with option -uniq.')
        # TODO: Change or document it
        # uniq is deprecated : to do a signle file, just give a trascript.pir
        # with only 1 sequence in it

    ###################### phylogenetic inference ######################

    if doPhylo:

        print("--------------------------------------------")
        print("Running phylogenetic reconstruction step...")
        print("--------------------------------------------")

        bestAssign = []
        bestScore = 10000

        # couts de changement d'etat des exons pour les especes ancestrales:
        # costs are ok
        # C1 : niveau proteique (etape 2)
        # C1= [ [0,cm,0], [cm,0,0], [0,cm,0] ]
        C1_1 = [[[0, cm, 0], [cm, 0, 0]], [[0, 0, cm], [0, 0, 0]]]
        # 1 -> 1                 0/2 -> 1
        C1_2 = [[0, 0, 0], [cm, 0, 0]]
        # 0->2     1/2->2
        C1_0 = [0, cm, 0]

        # C0 : niveau genique (etape 1)
        C0 = [[0, 1, 1], [1, 0, 1], [1, 1, 0]]

        priority = (1, 0, 2)

        costs = (CB, CD, cm)
        costMat = (C0, C1_0, C1_1, C1_2)

        dat, AllExons = initData.initTree(inputTree, inputFile, prune)
        print("The exons are:")
        print(AllExons)
        # nExons = range(1, len(AllExons)+1)
        # print nExons
        # nExons = [str(i) for i in nExons]
        # print(type(nExons))

        if topoStart == {}:
            res = ip.bestTopology(dat, AllExons, nbIt, costs, costMat,
                                  priority, SUFF, initBest, slowMode,
                                  topoStart, withMemory, outputDir)
        else:
            if printOnly:
                res = [topoStart]
            else:
                res = ip.bestWideTopology(dat, AllExons, costs, costMat,
                                          priority, SUFF, initBest, topoStart)

        if len(res[0].nodes()) > 0:
            exSt = ip.ex_state(ip.exState(dat, costMat, AllExons), priority)
            distTabs = ip.leafScoreTabs(dat, exSt, costMat, AllExons)
            if printOnly:
                tree = topoStart
                ip.writeOutput((tree, 0, 0), exSt, SUFF,costs, AllExons,
                outputDir)
            else:
                tree = ip.leafAssign(res[0], exSt, distTabs, costMat, AllExons)
                ip.writeOutput((tree, res[2], res[3]), exSt,SUFF, costs, AllExons,
                outputDir)
        else:
            print("No suitable topology could be found.")

    ###################### molecular modeling ######################
    if doModel:
        print("--------------------------------------------")
        print("Running molecular modeling step...")
        print("--------------------------------------------")
        HHBLITS, ADDSS, HHMAKE, HHSEARCH, HHMODEL, CONTEXTLIB = mi.getProgramPath(
        HHLIB)
        os.chdir(outputDir)

        # create as many directories as fasta input files
        # and inside as many fasta files as transcripts
        if not uniq and not onlyQuality:
            print("prepare intputs...")
            print(pathTransSeqs, outputDir)
            pathlib.Path(outputDir).mkdir(parents=True, exist_ok=True)
            for root, dirs, files in os.walk(pathTransSeqs, topdown=False):
                for name in dirs:
                    pathlib.Path(os.path.join(outputDir, name)).mkdir(parents=True, exist_ok=True)
                    mi.parseFromThorAxe(root+'/'+name, os.path.join(outputDir, name))
        # determine the number of templates that will be retained
        selTemp = ""
        for i in range(nbTemp):
            selTemp += str(i+1)
            selTemp += " "
        # NOTE: for nbTemp == 5 this block generates '1 2 3 4 5 '

        # perform the whole 3D modelling process

        # the 'if' is to perform the modelling in __pycache__
        # might add a list of forbidden names for less errors
        dirs = [x[0] for x in os.walk('.') if not '__pycache__' in x[0]]
        # NOTE: dirs is at leat ['.'] and it has subfolders,
        # e.g. ['.', './folder', './folder/subfolder']

        if not onlyQuality:
            print("launch the 3D modelling process")
            if uniq:
                if not os.path.isfile(pathTransSeqs):
                    raise ValueError('--instruct should be a file if '
                                     '--uniq is used')

                mydir, trans = os.path.split(pathTransSeqs)
                os.chdir(mydir)
                mi.runModelProcess(HHBLITS, ADDSS, HHMAKE, HHSEARCH, HHMODEL,
                                   HHDB, STRUCTDB, ALLPDB, NCPU, trans,
                                   selTemp, only3D, CONTEXTLIB)
                os.chdir(outputDir)
            else:
                print("Molecular modelling in: {}".format(str(dirs[1:])))
                for mydir in dirs[1:]:
                    os.chdir(mydir)
                    for trans in glob.glob('*.fa') + glob.glob('*.fasta') +\
                            glob.glob('*.faa'):
                        stop = False
                        name = os.path.splitext(trans)[0]
                        it = 1
                        while not stop:
                            mi.runModelProcess(HHBLITS, ADDSS, HHMAKE, HHSEARCH,
                                               HHMODEL, HHDB, STRUCTDB, ALLPDB,
                                               NCPU, trans, selTemp, only3D,
                                               CONTEXTLIB, it)
                            stop = mi.parseFromPirAnnotated(name, it)
                            it += 1
                            if it ==20: # stop with the infinite loop, temporary solution
                                stop = True
                        with open("number_of_loops.txt", "w") as f:
                            f.write(str(it-1))
                    os.chdir(outputDir)

        # assess the quality of the models
        # print('assess the quality of the models')
        # with open('quality.sum', 'w') as fOUT:
        #     fOUT.write('# procheck: Ideally, scores should be above -0.5. '
        #                'Values below -1.0 may need investigation.\n')
        #     fOUT.write('# dope: This is a Z-score; positive scores are '
        #                'likely to be poor models, while scores lower '
        #                'than -1 or so are likely to be native-like.\n')
        #     fOUT.write('transcript\tlenFull\tpercentSS\tlenModel\tdihedrals\t'
        #                'covalent\toverall\tdope\trSurf\trHydroph\n')
        #     for mydir in dirs[1:]:
        #         try:
        #             os.chdir(mydir)
        #             utils.makedirifnot('procheck')
        #             os.chdir('procheck')
        #             for prot in glob.glob('../*.B99990001.pdb'):
        #                 mydir, trans = os.path.split(prot)
        #                 pref = trans.split('.')[0]
        #                 path_pref = os.path.join('..', pref)
        #                 lenModel = mi.computeLenModel(path_pref)
        #                 lenFull, percentSS = mi.computeSS(path_pref)
        #                 dihedrals, covalent, overall = mi.assessQuality(
        #                     PROCHECK, prot, pref)
        #                 zscore = mi.assessNormalizedDopeScore(prot)
        #                 rSurf, rHydroph = mi.computeRatioSASA(
        #                     NACCESS, prot, pref)
        #                 fOUT.write(pref + '\t' + str(lenFull) + '\t' +
        #                            str(percentSS) + '\t' + str(lenModel) +
        #                            '\t' + dihedrals + '\t' + covalent +
        #                            '\t' + overall + '\t' + str(zscore) +
        #                            '\t' + str(rSurf) +
        #                            '\t' + str(rHydroph) + '\n')
        #             os.chdir('..')
        #             utils.clear_folder('procheck')
        #             os.rmdir('procheck')
        #             os.chdir(outputDir)
        #         except:
        #             print('Problem with ' + mydir)


def main():
    start = time.time()
    args = parse_command_line()
    doit(args.phylo,
         args.model,
         args.tree,
         args.transcripts,
         args.b,
         args.d,
         args.instruct,
         args.m,
         args.ni,
         args.nt,
         args.outputdir,
         args.only3D,
         args.onlyquality,
         args.s,
         args.suffix,
         args.topo,
         args.uniq,
         args.printonly,
         args.withmemory,
         not args.noprune,
         args.hhlib,
         args.hhdb,
         args.structdb,
         args.allpdb,
         args.ncpu
    )
    end = time.time()
    print("finished in {} seconds".format(end - start))

if (__name__ == '__main__'):
    main()
