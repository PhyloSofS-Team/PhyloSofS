# -*- coding: utf-8 -*-

# Copyright (c) 2014-2018: Elodie Laine, Hugues Richard, Adel Ait-hamlat
# and Diego Javier Zea.
# This code is part of the phylosofs package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

# Modele d'inference de phylogenies de transcrits
# Pour etudier l'apparition et la fixation
# d'evenements d'epissage au cours de l'evolution

from __future__ import division
import sys
import os
import glob
import random
import initData
import pickle as pk
import utils
import inferPhylo as ip
import modelIsoforms as mi
import argparse

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
            parser.error("phylosofs requires " + group +
                         " if " + argument + " is used.")
    else:
        if required and arg_dict[arg_name] is None:
            parser.error("phylosofs requires " + argument +
                         " if " + group + " is used.")
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
        If you use it, please cite:

        Ait-hamlat A, Polit L, Richard H, Laine E. Transcripts evolutionary
        conservation and structural dynamics give insights into the role of
        alternative splicing for the JNK family. bioRxiv. 2017 Jan 1:119891.

        It has been developed at LCQB (Laboratory of Computational and
        Quantitative Biology), UMR 7238 CNRS, Sorbonne Université.
        """,
    )

    parser.add_argument(
        '-P', '--phylo',
        help='Do the phylogenetic inference',
        action='store_true'
    )
    parser.add_argument(
        '-M', '--model',
        help='Do the molecular modelling',
        action='store_true'
    )
    parser.add_argument(
        '--inseq',
        help='text file containing input data: string representing '
        'the gene tree in Newick format on the first line, and then '
        'the list of transcripts for each leave (leaf_name: t1,t2,t3...)'
    )
    parser.add_argument(
        '-c',
        help='text file containing values for parameters '
        '(for molecular modeling part)'
    )  # REQUIRED
    parser.add_argument(
        '-b',
        help='birth cost',
        default=5
    )
    parser.add_argument(
        '-d',
        help='death cost',
        default=3
    )
    parser.add_argument(
        '--instruct',
        help='text file containing input data: either a directory where '
        'multifasta files are located (one file per species) or a single '
        'fasta file with the sequence of only one transcript (in that '
        'case, the unique option must be set to TRUE)',
        default='.'
    )
    parser.add_argument(
        '-m',
        help='mutation cost',
        default=2
    )
    parser.add_argument(
        '--ni',
        help='number of iterations',
        default=1
    )
    parser.add_argument(
        '--nt',
        help='number of templates to retain',
        default=5
    )
    parser.add_argument(
        '--noprune',
        help='disable the removal of exons that appear in only one transcript',
        action='store_true'
    )
    parser.add_argument(
        '-o',
        help='output directory',
        default='.'
    )
    parser.add_argument(
        '--only3D',
        help='perform only the 3D modeling step (skip template search)',
        action='store_true'
    )
    parser.add_argument(
        '--onlyquality',
        help='perform only the 3D models quality assessment',
        action='store_true'
    )
    parser.add_argument(
        '--printonly',
        help='perform only the generation of the PDF file enabling to '
        'visualize the transcripts',
        action='store_true'
    )
    parser.add_argument(
        '--uniq',
        help='treat only one transcript whose sequence is taken from '
        'the fasta file indicated by the -i option',
        action='store_true'
    )
    parser.add_argument(
        '-s',
        help='starting score (by default: not considered), '
        'if no score is given, the algorithm starts the search '
        'by the topology corresponding to the maximum number of '
        'binary subnodes at each nodes (forest with the smallest '
        'possible number of trees), otherwise it starts from a '
        'randomly chosen topology'
    )
    parser.add_argument(
        '--suffix',
        help='suffix, it is _(birth)(death)(mutation)_(iterations) by '
        'default, e.g. _532_1'
    )
    parser.add_argument(
        '--topo',
        help="initial topology (by default: maximum or random topology), "
        "or transcripts' phylogeny to be printed out "
        "(if the -printonly option is active)"
    )
    parser.add_argument(
        '--withmemory',
        action='store_true'
    )

    args = parser.parse_args()

    if args.phylo is None and args.model is None:
        parser.error("phylosofs requires --phylo or/and --model flags.")

    arg_dict = vars(args)

    check_argument_groups(parser, arg_dict, '--phylo', '--inseq', True)
    check_argument_groups(parser, arg_dict, '--model', '-c', True)

    # Check flag arguments
    if not args.model:
        if args.only3D:
            parser.error("phylosofs requires --model if --only3D is used.")

    return args


def main(doPhylo,
         doModel,
         inputFile,
         configFile,
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
         prune
         ):
    """
    main(args.phylo,
         args.model,
         args.inseq,
         args.c,
         args.b,
         args.d,
         args.instruct,
         args.m,
         args.ni,
         args.nt,
         args.o,
         args.only3D,
         args.onlyquality,
         args.s,
         args.suffix,
         args.topo,
         args.uniq,
         args.printonly,
         args.withmemory,
         not args.noprune
         )
    """
    random.seed()

    if doPhylo:
        if not os.path.isfile(inputFile):
            sys.stderr.write("You must give a valid input file for "
                             "phylogenetic inference.")

    if doModel:
        if os.path.isfile(configFile):
            HHBLITS, ADDSS, HHMAKE, HHSEARCH, HHMODEL, PROCHECK, NACCESS,\
                HHDB, STRUCTDB, ALLPDB, NCPU, CONTEXTLIB = mi.init(configFile)
        else:
            sys.stderr.write("You must give an existing configuration file "
                             " for molecular modeling. See usage instructions")

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

        dat, AllExons = initData.initTree(inputFile, prune)
        print "The exons are:"
        print AllExons
        nExons = range(len(AllExons))
        print nExons

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
                ip.writeOutput((tree, 0, 0), exSt, SUFF,
                               costs, AllExons, outputDir)
            else:
                tree = ip.leafAssign(res[0], exSt, distTabs, costMat, AllExons)
                ip.writeOutput((tree, res[2], res[3]), exSt,
                               SUFF, costs, AllExons, outputDir)
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
            mi.prepareInputs(pathTransSeqs)

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
                mi.runModelProcess(HHBLITS, ADDSS, HHMAKE, HHSEARCH,
                                   HHMODEL, HHDB, STRUCTDB, ALLPDB, NCPU,
                                   # './'+trans, 
                                   selTemp, only3D)
                os.chdir('..')
            else:
                print dirs
                for mydir in dirs[1:]:
                    os.chdir(mydir)
                    for trans in glob.glob('./*.fa'):
                        mi.runModelProcess(HHBLITS, ADDSS, HHMAKE, HHSEARCH,
                                           HHMODEL, HHDB, STRUCTDB, ALLPDB,
                                           NCPU, trans, selTemp, only3D)
                    os.chdir('..')

        # assess the quality of the models
        print 'assess the quality of the models'
        with open('quality.sum', 'w') as fOUT:
            fOUT.write(
                '# procheck: Ideally, scores should be above -0.5. '
                'Values below -1.0 may need investigation.\n')
            fOUT.write('# dope: This is a Z-score; positive scores are '
                       'likely to be poor models, while scores lower '
                       'than -1 or so are likely to be native-like.\n')
            fOUT.write(
                'transcript\tlenFull\tpercentSS\tlenModel\tdihedrals\t'
                'covalent\toverall\tdope\trSurf\trHydroph\n')
            for mydir in dirs[1:]:
                try:
                    os.chdir(mydir)
                    utils.makedirifnot('procheck')
                    os.chdir('procheck')
                    for prot in glob.glob('../*.B99990001.pdb'):
                        mydir, trans = os.path.split(prot)
                        pref = trans.split('.')[0]
                        lenModel = mi.computeLenModel('../'+pref) # TODO : Test on windows
                        lenFull, percentSS = mi.computeSS('../'+pref)
                        dihedrals, covalent, overall = mi.assessQuality(
                            PROCHECK, prot, pref)
                        zscore = mi.assessNormalizedDopeScore(prot)
                        rSurf, rHydroph = mi.computeRatioSASA(
                            NACCESS, prot, pref)
                        fOUT.write(pref + '\t' + str(lenFull) + '\t' +
                                   str(percentSS) + '\t' + str(lenModel) +
                                   '\t' + dihedrals + '\t' + covalent +
                                   '\t' + overall + '\t' + str(zscore) +
                                   '\t' + str(rSurf) + 
                                   '\t' + str(rHydroph) + '\n')
                    os.chdir('..')
                    utils.clear_folder('procheck')
                    os.rmdir('procheck')
                    os.chdir('..')
                except:
                    print 'Problem with ' + mydir


if (__name__ == '__main__'):

    args = parse_command_line()

    main(args.phylo,
         args.model,
         args.inseq,
         args.c,
         args.b,
         args.d,
         args.instruct,
         args.m,
         args.ni,
         args.nt,
         args.o,
         args.only3D,
         args.onlyquality,
         args.s,
         args.suffix,
         args.topo,
         args.uniq,
         args.printonly,
         args.withmemory,
         not args.noprune
         )
