# -*- coding: utf-8 -*-

# Copyright (c) 2014-2018: Adel Ait-hamlat, Elodie Laine, Lélia Polit
# and Diego Javier Zea.
# This code is part of the phylosofs package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

# module to reconstruct 3D models of isoforms

import glob
import os
import sys
import pathlib
import shutil
import subprocess
import warnings
import gzip
import pdb
import pkg_resources
import configparser
import numpy as np
from collections import OrderedDict
from Bio.PDB import *
_config = configparser.ConfigParser()

try:
    import modeller
    from modeller import automodel
    _MODELLER_MESSAGE = ""
except ImportError:
    warnings.warn('modeller is not installed', ImportWarning)
    _MODELLER_MESSAGE = "Please install modeller: " \
                        "https://salilab.org/modeller/"
except Exception as err:
    err_message = str(err)
    if err_message.startswith("check_lice_E"):
        warnings.warn('invalid modeller license key', ImportWarning)
        _MODELLER_MESSAGE = '\n'.join(err_message.split('\n')[1:])
    else:
        raise err

# ------------------ UPLOAD INIT VALUES ----------------------------- #
#
#


def init(configFile):

    _config.read(os.path.expanduser(configFile))
    # upload paths for librairies
    HHBLITS = _config.get("PROGRAMS", "HHBLITS")
    ADDSS = _config.get("PROGRAMS", "ADDSS")
    HHMAKE = _config.get("PROGRAMS", "HHMAKE")
    HHSEARCH = _config.get("PROGRAMS", "HHSEARCH")
    HHMODEL = _config.get("PROGRAMS", "HHMODEL")
    PROCHECK = _config.get("PROGRAMS", "PROCHECK")
    NACCESS = _config.get("PROGRAMS", "NACCESS")
    HHDB = _config.get("DATABASES", "HHDB")
    STRUCTDB = _config.get("DATABASES", "STRUCTDB")
    ALLPDB = _config.get("DATABASES", "ALLPDB")
    NCPU = _config.get("OPTIONS", "NCPU")
    CONTEXTLIB = _config.get("DATA", "CONTEXTLIB")

    return (HHBLITS, ADDSS, HHMAKE, HHSEARCH, HHMODEL, PROCHECK, NACCESS, HHDB,
            STRUCTDB, ALLPDB, NCPU, CONTEXTLIB)


# This function takes as an argument a directory path to the hh-suite
# and gets every executable paths


def getProgramPath(HHLIB):
    # HHBLITS, ADDSS, HHMAKE, HHSEARCH, HHMODEL, CONTEXTLIB = ''
    HHBLITS, HHMAKE, HHSEARCH, HHMODEL, CONTEXTLIB = ('', '', '', '', '')
    for root, dirs, files in os.walk(HHLIB, topdown=False):
        for name in files:
            program = os.path.join(root, name)
            if (program.split("/")[-2:] == ['bin', 'hhblits']):
                HHBLITS = program
            elif (program.split("/")[-2:] == ['bin', 'hhmake']):
                HHMAKE = program
            elif (program.split("/")[-2:] == ['bin', 'hhsearch']):
                HHSEARCH = program
            elif (program.split("/")[-1] == 'hhmakemodel.py'
                  and program.split("/")[-3] == 'build'):
                HHMODEL = program
            elif (program.split("/")[-1] == 'context_data.lib'):
                CONTEXTLIB = program
    if ((HHBLITS == '') or (HHMAKE == '') or (HHSEARCH == '')
            or (HHMODEL == '') or (CONTEXTLIB == '')):
        print(
            'Could not locate the HHSUITE directory, please check your --hhlib argument'
        )
        sys.exit(1)
    return (HHBLITS, HHMAKE, HHSEARCH, HHMODEL, CONTEXTLIB)


# sys.path.append(STRUCTURE)
# sys.path.append(SEQUENCE)
#
# -------------------------------------------------------------------- #

# read a multi fasta file
# if it is the initial file, a special treatment is applied to names


def readFastaMul(fic, init=False):

    res = []
    seq = ''
    name = ''
    go = False

    path, name = os.path.split(fic)
    ext = name.split('.')[1]

    with open(fic, 'r') as f:
        for line in f.readlines():
            if line.strip() != '':
                if line[0] == '>':
                    if seq != '':
                        res.append((name, seq))
                        seq = ''
                    if init:
                        if 'peptide' in line:
                            name = line[0] + line.split(":")[1]
                            name = name.split(" ")[0]
                            name = name.replace('>', '')
                            name = name.split(" ")[0] + "\n"
                            name = name.strip()
                            go = True
                        else:
                            go = False
                    else:
                        name = line.replace(">", "").replace("\n", "")
                        go = True
                else:
                    if ext == 'pir' and seq == '':
                        name += '\n' + line.strip()
                        seq = '?'
                    else:
                        if go:
                            if seq == '?':
                                seq = line.strip()
                            else:
                                seq = seq + line.strip()
        # add the last one
        if go:
            res.append((name, seq))
    return res


def parse_pir(file_path):
    """
    Parse an annotated pir file returning three dictionaries.

    The dictionary contains the sequences, their lenghts and the exon sequences.
    """
    sequences = {}
    lengths = {}
    exons_seqs = {}
    with open(file_path, encoding='utf-8') as pir:
        line_number = 0
        for line in pir:
            line_number += 1
            if line_number == 1 and not line.startswith(">"):
                raise Exception(
                    "{} is not a correct pir file, it should start with '>'".
                    format(file_path))
            if line.startswith(">"):
                header = line
            elif line.endswith("*\n"):
                sequences[header] = line
            else:
                # we can also use the 3rd column of the header, but this seems to be easier
                # and is fast
                exons_seqs[header] = line
                line = line.rstrip('\r\n')
                exons_lengths = []
                character = ""
                for i, char in enumerate(line):
                    if character != char:
                        exons_lengths.append(i + 1)
                        character = char
                # at the the lengths the end of the sequence
                exons_lengths.append(len(line))
                lengths[header] = exons_lengths
    return sequences, lengths, exons_seqs


# parse a 'transcripts.pir' file form the ThoxAxe output
# and transforms it to be usable easily by PhyloSofS
def parseFromThorAxe(pathTransSeqs, outputDir):
    (sequences, lengths,
     exons_seqs) = parse_pir(os.path.join(pathTransSeqs, "transcripts.pir"))
    for i in sequences:
        if "ENSG0" in i:  # take only the human transcripts
            j = '_'.join(i.split()[0:2]).replace('>P1;', '')
            pathlib.Path(outputDir + '/' + j).mkdir(parents=True,
                                                    exist_ok=True)

            with open(outputDir + '/' + j + '/' + j + '.fasta', 'w+') as file:
                file.write(i)
                # remove the "*\n" at the end
                file.write(sequences[i][:-2] + '\n')

            with open(outputDir + '/' + j + '/' + j + '_1.fa', 'w+') as file:
                file.write(i)
                # remove the "*\n" at the end
                file.write(sequences[i][:-2] + '\n')

            with open(outputDir + '/' + j + '/' + j + '.exons_lengths.txt',
                      'w+') as file:
                file.write(i)
                file.write(' '.join(str(j) for j in lengths[i]) + '\n')

            with open(outputDir + '/' + j + '/' + j + '_1.exons_lengths.txt',
                      'w+') as file:
                file.write(i)
                file.write(' '.join(str(j) for j in lengths[i]) + '\n')

            with open(outputDir + '/' + j + '/' + j + '_annotated.pir',
                      'w+') as file:
                file.write(i)
                file.write(exons_seqs[i])
                file.write(sequences[i])
                file.write(len(exons_seqs[i][:-1]) * '0')


# Parse the sequence from the output pir with the annotations on what to take
# for the iterative version of hhblits:
# 0: never ran, so no templates found
# 1: hhblits ran, templates found
# 2: hhblits ran once, no templates found (or bad ones)
# 3: hhblits ran a second time with another fragment, still no templates found --> stop iteration
# The function returns True if the iterations should stop, otherwise returns False.


def parseFromPirAnnotated(name, it):
    go = False
    while go == False:
        go = True
        with open("./" + name + "_annotated.pir", "r") as f:
            lines = f.readlines()
        try:
            first = next(index for index, value in enumerate(lines[3])
                         if (value == '0' or value == '2'))
        except:  # no 0 on 2 in the sequence, no more runs to do
            return True
        for i in range(first, len(lines[3])):
            if all([lines[3][i] != '0', lines[3][i] != '2']):
                last_contiguous_0 = i
                break
            last_contiguous_0 = len(lines[3])
        a = np.array(list(lines[1][first:last_contiguous_0]))
        header = ' '.join(lines[0].split()[0:2]) + ' ' + \
            ''.join(OrderedDict.fromkeys(a).keys())

        exons_lengths = []
        character = ""
        line = lines[1][first:last_contiguous_0]
        for i in range(len(line)):
            if character != line[i]:
                exons_lengths.append(i + 1)
                character = line[i]
        # at the the lengths the end of the sequence
        exons_lengths.append(len(line) + 1)

        if (len(line) < 5):
            # modify the pir_annotated
            with open(name + '_annotated.pir', 'w') as f:
                f.write(lines[0])
                f.write(lines[1])
                f.write(lines[2])
                f.write(lines[3][0:first] + "3" + lines[3][last_contiguous_0:])
            go = False

        if name + '.fasta' in os.listdir():
            os.remove(name + '.fasta')
        elif name + '.exons_lengths.txt' in os.listdir():
            os.remove(name + '.exons_lengths.txt')

        # writing of the fasta
        with open(name + '.fasta', 'w') as f:
            f.write(header + '\n')
            f.write(lines[2][first:last_contiguous_0])

        with open(name + '_' + str(it + 1) + '.fa', 'w') as f:
            f.write(header + '\n')
            f.write(lines[2][first:last_contiguous_0])

        # writing of the exons_lengths
        with open(name + '.exons_lengths.txt', 'w') as f:
            f.write(header + '\n')
            f.write(' '.join(str(j) for j in exons_lengths) + '\n')

        with open(name + '_' + str(it + 1) + '.exons_lengths.txt', 'w') as f:
            f.write(header + '\n')
            f.write(' '.join(str(j) for j in exons_lengths) + '\n')
        print(go)
        if go != False:
            if len(lines[2][first:last_contiguous_0]) < 5:
                return True
            else:
                return False


# a partir d'une banque de transcrit, ecrit chaque transcrit dans un fichier
# en ne retenant que son nom et la sequence

# def parse(fic):
#     seqs = readFastaMul(fic, True)
#     # nameDir = fic.split("/")
#     # nameDir = nameDir[len(nameDir)-1].split(".")[0]
#     nameDir = os.path.basename(os.path.dirname(fic)).split(".")[0]
#     if not os.path.exists(nameDir):
#         os.makedirs(nameDir)
#     os.chdir(nameDir)
#     writeFastas(seqs)
#     fics = glob.glob('*.fa')
#     os.chdir("..")
#     return fics


def writeFastas(seqs):
    for seq in seqs:
        # there is many replaces to avoid problems later, with the JPred API for example (allows only letters, numbers and underscores in the name)
        with open((seq[0].replace('/', '_').replace(" ", "_").replace(
                "|", "_").replace("-", "_").replace("=", "_"))[0:40] + ".fa",
                  'w') as f:
            f.write('>' + seq[0] + '\n')
            i = 0
            while i < len(seq[1]):
                f.write(seq[1][i:(i + 60)])
                f.write('\n')
                i += 60
    return 1


def treatAli(fic):
    seqs = readFastaMul(fic)
    leseq = len(seqs[0][1])
    n = len(seqs[0][1])
    nbSeqs = len(seqs)
    borders = [0, 0]
    valInit = [-1, n - 1]
    valInc = [1, -1]
    for k in range(2):
        i = valInit[k]
        allTirets = True
        while allTirets:
            i += valInc[k]
            allTirets = True
            s = 0  # NOTE: This was 1, but that created an infinite loop
            # when nbSeqs was 1... Was there a reson for s = 1?
            while s < nbSeqs and allTirets:
                allTirets = (seqs[s][1][i] == '-')
                s += 1
        borders[k] = i
    res = writePirMul(fic, seqs, borders)
    return res


def writePirMul(fic, seqs, borders=[]):
    if borders == []:
        borders = [0, len(seqs[0][1])]
    else:
        secondLine = seqs[0][0].split('\n')[1].split(':')
        start = int(secondLine[2]) + borders[0]
        startstr = "% 4d" % start
        end = int(secondLine[4]) - (len(seqs[0][1]) - borders[1]) + 2
        endstr = "% 4d" % end
        newNam = seqs[0][0].split('\n')[0] + '\n' + secondLine[0] + ':' + \
            secondLine[1] + ':' + startstr + ':' + secondLine[3] + ':' + endstr
        for j in range(5, len(secondLine)):
            newNam += ':' + secondLine[j]
        newSeq = (newNam, seqs[0][1])
        seqs[0] = newSeq

    nam = fic.split('.')[0]
    # os.system('mv ' + fic + ' ' + fic + '_old')
    shutil.move(fic, fic + '_old')
    with open(fic, 'w') as f:
        for seq in seqs:
            f.write('>' + seq[0] + '\n')
            i = 0
            s = seq[1][borders[0]:(borders[1] + 1)] + '*'
            while i < len(s):
                f.write(s[i:(i + 100)])
                f.write('\n')
                i += 100
    return [int(secondLine[2]), int(secondLine[4]), start, end]


# extract each transcript sequence

# def prepareInputs(pathTransSeqs):
#     files = glob.glob('*.fa') + glob.glob('*.fasta') + glob.glob('*.faa')
#     for fic in files:
#         # nameDir = fic.split("/")
#         # nameDir = nameDir[len(nameDir)-1].split(".")[0]
#         nameDir = os.path.basename(os.path.dirname(fic)).split(".")[0]
#         ficsOUT = parse(fic)
#         for f in ficsOUT:
#             nameRTF = os.path.join(pathTransSeqs, f.split(".")[0]+".rtf")
#             # nameRTF = '.'.join(f.split('.')[:-1]) + ".rtf"
#             if os.path.isfile(nameRTF):
#                 # os.system("cp "+nameRTF+" "+nameDir+"/")
#                 shutil.copy2(nameRTF, nameDir)
#     os.chdir(folder)


def prepareInputs(transcriptsdir, outputdir):
    here = os.getcwd()  # To be able to go back...
    if os.path.isdir(transcriptsdir):
        os.chdir(transcriptsdir)
        files = glob.glob('*.fa') + glob.glob('*.fasta') + glob.glob('*.faa')
        os.chdir(outputdir)
        for file in files:
            path = os.path.join(transcriptsdir, file)
            seqs = readFastaMul(path, False)
            namedir = file.split('.')[0]
            if not os.path.exists(namedir):
                os.mkdir(namedir)
            os.chdir(namedir)
            writeFastas(seqs)
            ffiles = glob.glob('*.fa')
            for ff in ffiles:
                namertf = os.path.join(transcriptsdir,
                                       ff.replace('.fa', '') + ".rtf")
                if os.path.isfile(namertf):
                    shutil.copy2(namertf, os.path.join(outputdir, namedir))

            os.chdir(outputdir)
    else:
        print(transcriptsdir + " is not a folder (running in " + here + ").")
    os.chdir(here)


# ecrit dans un fichier
def writeToFile(name, chaine):
    with open(name, 'w') as n:
        n.write(chaine)
    return 0


# read a pdb file (IN) and split the different chains in several pdb files


def splitChainsPDB(fic, code, ext):
    rfile = open(fic, "r")

    # os.system("rm "+fic+"_* 2> err.txt")
    for file in glob.glob(fic + "_*"):
        os.remove(file)

    lines = rfile.readlines()
    nbLines = len(lines)

    ligne = ""
    count = 0

    while count < nbLines:

        while not ligne.startswith("ATOM") and (count < nbLines):
            try:
                ligne = lines[count]
            except:  # TODO : Try not to use try here or use a better except.
                print("can't read line")
            count += 1

        if count < nbLines:
            chain = ligne[21]
            ficOUT = code + "_" + chain + '.' + ext
            wfile2 = open(ficOUT, "w")
            wfile2.writelines(ligne)
            try:
                ligne = lines[count]
            except:  # TODO : Try not to use try here or use a better except.
                print("can't read line")
            count += 1
            while (ligne[21] == chain) and ((ligne.startswith("ATOM")) or
                                            (ligne.startswith("HETATM")) or
                                            (ligne.startswith("ANISOU"))) and \
                  (count < nbLines):
                wfile2.writelines(ligne)
                try:
                    ligne = lines[count]
                except:  # TODO : Try not to use try or use a better except.
                    print("can't read line")
                count += 1
            wfile2.close()


def assessNormalizedDopeScore(pdb):
    if _MODELLER_MESSAGE != "":
        raise ImportError(_MODELLER_MESSAGE)

    env = modeller.environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')  # TODO: Test on Windows
    env.libs.parameters.read(file='$(LIB)/par.lib')
    # Read a model previously generated by Modeller's automodel class
    mdl = modeller.scripts.complete_pdb(env, pdb)
    zscore = mdl.assess_normalized_dope()
    return zscore


def assessQuality(PROCHECK, prot, pref):
    # os.system(PROCHECK + ' ' + prot + ' 1.5')
    subprocess.call([PROCHECK, prot, ''])
    with open(pref + '.B99990001.sum', 'r') as f:
        lines = f.readlines()
        for line in lines:
            fields = line.split()
            if 'G-factors' in fields:
                dihedrals, covalent, overall = fields[3], fields[5], fields[7]
    # os.system('mv ' + pref + '.B99990001.sum ../')
    shutil.move(pref + '.B99990001.sum', '..')
    return dihedrals, covalent, overall


def computeSS(pref):
    seqs = readFastaMul(pref + '.a3m')
    lenFull = len(seqs[0][1])
    count = 0
    for i in range(lenFull):
        if seqs[0][1][i] != 'C' and int(seqs[1][1][i]) > 5:
            count += 1
    return lenFull, float(count) / lenFull


def computeLenModel(pref):
    seqs = readFastaMul(pref + '.pir')
    count = 0
    for let in seqs[0][1]:
        if let != '-' and let != '*':
            count += 1
    return count


def computeRatioSASA(NACCESS, prot, pref):
    # os.system('cp ' + prot + ' tmp.pdb')
    shutil.copy2(prot, 'tmp.pdb')
    # os.system(NACCESS + ' tmp.pdb')
    subprocess.call([NACCESS, 'tmp.pdb'])
    hydroph = ["VAL", "ILE", "LEU", "MET", "PHE", "TRP", "CYS"]
    with open('tmp.rsa', 'r') as f:
        lines = f.readlines()
        n = 0
        nsurf = 0
        nhydophsurf = 0
        nhydroph = 0
        for line in lines:
            if line[0:3] == "RES":
                fields = line.split()
                n += 1
                restype = fields[1]
                rsa = fields[4]
                if float(rsa) > 25:
                    nsurf += 1.0
                    if restype in hydroph:
                        nhydophsurf += 1.0
                        nhydroph += 1
                else:
                    if restype in hydroph:
                        nhydroph += 1
    # os.system('mv tmp.rsa ../' + pref + '.rsa')
    shutil.move('tmp.rsa', os.path.join('..', pref + '.rsa'))
    # os.system('rm tmp.*')
    for file in glob.glob('tmp.*'):
        os.remove(file)

    return nsurf / n, nhydophsurf / nhydroph


def model3D(fic, ALLPDB, pdb_extension='.cif'):
    if _MODELLER_MESSAGE != "":
        raise ImportError(_MODELLER_MESSAGE)

    seqs = readFastaMul(fic)

    if len(seqs) < 2:
        raise Exception("There aren't template sequences in %s." % fic)

    seq = seqs[0][0].split('\n')[0]
    seq = seq.split(';')[1]

    modeller.log.verbose()  # request verbose output
    env = modeller.environ()  # create a new MODELLER environment to build ...
    # ... this model in
    env.io.atom_files_directory = ['.']  # ['.', ALLPDB]
    knowns = []
    for i in range(1, len(seqs)):
        tmp = seqs[i][0].split('\n')[0]
        tmp = tmp.split(';')[1]
        knowns.append(tmp)
        kn = tmp.split('_')[0].lower()  # upper ?
        base_name = kn + pdb_extension
        nam = base_name  # + '.gz'
        if not os.path.isfile(kn + pdb_extension):
            # shutil.copyfile(ALLPDB + nam, "./" + nam)
            shutil.copy2(os.path.join(ALLPDB, nam), nam)
            # # os.system("gunzip ./" + nam)
            # with gzip.open(nam, 'rb') as f_in:
            #     with open(base_name, 'wb') as f_out:
            #         shutil.copyfileobj(f_in, f_out)

            # splitChainsPDB('pdb' + kn + pdb_extension, kn, 'pdb')
    knowns = tuple(knowns)
    a = automodel.automodel(env, alnfile=fic, knowns=knowns, sequence=seq)
    a.max_molpdf = 1e12
    a.starting_model = 1  # index of the first model
    a.ending_model = 1  # index of the last model
    # (determines how many models to calculate)
    cdir = os.getcwd()
    print("cdir = ", cdir)
    a.make()  # do the actual homology modeling

    return 1


def colorBFactor(name):

    with open("./" + name + "_annotated.pir") as f:
        lines = f.readlines()
    exons = lines[1]
    char = ""
    ex = []
    exon_color = []
    liste = [1, 4, 2, 5, 3, 7]

    for i in range(len(exons)):
        if exons[i] != char:
            char = exons[i]
            ex.append(char)
    for i in range(len(exons)):
        index = int(ex.index(exons[i]))
        exon_color.append(liste[(index) % len(liste)])

    parser = PDBParser()
    structure = parser.get_structure('str', name + ".B99990001.pdb")
    for residue in structure.get_residues():
        for atom in residue:
            atom.bfactor = exon_color[residue.id[1] - 1]

    io = PDBIO()
    io.set_structure(structure)
    io.save(name + '_colored.pdb')


def readRTF(filename):
    f = open(filename, "rb")
    lines = f.readlines()
    f.close()
    s = ""
    for line in lines:
        s = s + line[:-1]

    Exons = []
    count = 0
    currCol = ""
    tmp = s.split(r"{\cf")
    for i in tmp[1:]:
        i = i.replace("\\\'46", "F")
        col = i[0]
        if col != currCol:
            currCol = col
            if col != "7":
                count = count + 1
        end = i.find("}")
        if col != "7":
            Exons = Exons + ([count] * len(i[1:end]))
        else:
            Exons = Exons + ([-1] * len(i[1:end]))
    return Exons


# insert information of the exons on the occupancy column


def insertExons(myExons, trans):

    with open(trans + ".B99990001.pdb", 'r') as fIN:
        lines = fIN.readlines()
        with open(trans + "_ex.pdb", "w") as fOUT:
            count = -1
            currRes = ""
            for line in lines:
                if line[:4] == "ATOM":
                    residue = line[17:26]
                    if currRes != residue:
                        currRes = residue
                        count = count + 1
                    fOUT.write(line[:54] + "%6.2f" % myExons[count] +
                               line[60:])
                else:
                    fOUT.write(line)

    return 1


# annotate the generated PDBs with exon information


def annotate(trans, borders):
    try:
        # tmp = trans[2:].split('.')[0]
        tmp = os.path.splitext(trans)[0]
        rtf_file = tmp + '.rtf'
        if os.path.isfile(rtf_file):
            myExons = readRTF(tmp + '.rtf')
            start = borders[2] - borders[0]
            end = len(myExons) - (borders[1] - borders[3] + 1)
            myExons = myExons[start:end]
            insertExons(myExons, tmp)
        else:
            warnings.warn('There is not ' + rtf_file + ' file!')

        res = 1
    except Exception as e:
        print('Error: could not annotate the 3D model for ' + tmp + ':\n    ' +
              str(e))
        res = 0
    return res


def modify_query_name(pirfile):
    # changes the name "UKNP" of the query in the .pir file
    # hhmakemodel.py puts 'UKNP as the sequence name if it does not
    # find a structure file for it, resulting on a 'UNKP.pdb' by modeller
    # Here we change the name beforehand to avoid this issue
    with open(pirfile, 'r') as f:
        filedata = f.read()
    newname = pathlib.Path(pirfile).stem
    newdata = filedata.replace("UKNP", newname)
    with open(pirfile, 'w') as f:
        f.write(newdata)


def run_external_program(command_list):
    """
    Use subprocess to print and run the command_list.

    >>> run_external_program(["python",  "-c", "None"])
    <BLANKLINE>
    ┌───────────────────────────────────────────────────────────────────────────
    │ RUNNING: python -c None
    <BLANKLINE>
    <BLANKLINE>
    │ FINISHED WITH EXIT CODE: 0
    └───────────────────────────────────────────────────────────────────────────
    <BLANKLINE>
    0
    >>> run_external_program(["python",  "-c", "print(example)"])
    <BLANKLINE>
    ┌───────────────────────────────────────────────────────────────────────────
    │ RUNNING: python -c print(example)
    ...
    │ FINISHED WITH EXIT CODE: 1
    └───────────────────────────────────────────────────────────────────────────
    <BLANKLINE>
    1
    """
    str_command = subprocess.list2cmdline(command_list)
    horizontal_bar = "─" * 75
    print("")
    print("┌%s" % (horizontal_bar))
    print("│ RUNNING: %s" % (str_command))
    print("")
    exit_code = subprocess.call(command_list)
    print("")
    print("│ FINISHED WITH EXIT CODE: %s" % (exit_code))
    print("└%s" % (horizontal_bar))
    print("")
    return exit_code


_JULIA_SCRIPT = pkg_resources.resource_filename('phylosofs',
                                                'src/plots_with_exons.jl')
_GET_PDBS = pkg_resources.resource_filename('phylosofs', 'src/get_pdbs.jl')


def runModelProcess(HHBLITS, HHMAKE, HHSEARCH, HHMODEL, HHDB, STRUCTDB, ALLPDB,
                    NCPU, trans, selTemp, only3D, CONTEXTLIB, it, JULIA):
    # try:
    # tmp = trans[2:].split('.')[0]  # 'example.fa' --> 'ample'
    tmp = os.path.splitext(trans)[0]  # 'example.fa' --> 'example'
    if not only3D:
        pass
        # look for homologous and create the MSA alignment
        run_external_program([
            HHBLITS,
            "-cpu",
            NCPU,
            "-i",
            trans,
            "-d",
            HHDB,
            "-neffmax",
            "11",  # skip further search iterations when diversity Neff of query MSA becomes larger than neffmax
            "-all",  # show all sequences in result MSA; do not filter result MSA
            "-id",
            "100",  # maximum pairwise sequence identity
            "-cov",
            "80",  # minimum coverage with master sequence (%)
            "-oa3m",
            tmp + ".a3m",
            "-n",
            "3"
            # "-maxfilt", "50000"
        ])

        # generate a hidden Markov model (HMM) from the MSA
        run_external_program([
            HHMAKE,
            "-i",
            tmp + ".a3m",
            # minimum coverage with master sequence (%)
            "-cov",
            "80",
            "-id",
            "100",  # maximum pairwise sequence identity
            "-neff",
            "11"
        ])

        # search for homologs
        run_external_program([
            HHSEARCH,
            "-cpu",
            NCPU,
            "-v",
            "1",
            "-i",
            tmp + ".hhm",
            "-d",
            STRUCTDB,
            "-o",
            tmp + "_" + str(it) + ".hhr",
            "-p",
            "20",
            "-Z",
            "100",
            "-B",
            "100",
            "-seq",
            "1",
            "-aliw",
            "80",
            "-loc",
            "-maxres",
            "37000",
            "-ssm",
            "2",
            "-norealign",
            "-sc",
            "1",
            "-cov",
            "80",  # minimum coverage with master sequence (%)
            "-all",  # show all sequences in result MSA; do not filter result MSA
            "-id",
            "100",  # maximum pairwise sequence identity
            "-cs",
            CONTEXTLIB
        ])

        run_external_program([
            JULIA, _GET_PDBS, "--input", tmp + "_" + str(it) + ".hhr", "--pdb",
            ALLPDB
        ])

        # create the alignment for MODELLER and change the query name
        run_external_program([
            "python3",
            HHMODEL,
            "-noduplicates",
            tmp + "_" + str(it) + ".hhr",
            ALLPDB,
            tmp + "_" + str(it) + ".pir",
            "./",
            # "-m", "1","2","3","4","5","6","7","8","9","10"
        ])

        # modify the query name in the alignment file
        # 'UKNP' --> tmp_name
        modify_query_name(tmp + "_" + str(it) + ".pir")

    # split the alignment for a template summary
    # for i in range(1,6):
    #     run_external_program(["python3",
    #             HHMODEL, tmp + ".hhr",
    #             ALLPDB,
    #             tmp + "_" + str(i) + ".pir",
    #             "./",
    #             "-m", str(i)
    #             ])

    # treat the alignment file to remove N- and C-terminal loops
    # borders = treatAli(tmp + '.pir')

    # analysis of the templates
    run_external_program([JULIA, "--inline=no", _JULIA_SCRIPT, tmp, str(it)])

    # Create files for secondary structures and solvent accessibility using JPred 4 API
    # run_external_program(["python",
    # "../../jpred_api_test.py", tmp+".pir"])

    # generate the 3D models with Modeller
    try:
        # model3D(tmp + "_"  + str(it) + '.pir', ALLPDB)
        # annotate(trans, borders)
        res = 0
    except Exception as e:
        print('Error: could not build the 3D model for ' + tmp)
        print("reason : {}".format(e))
        res = 1
    # res = 0

    return res
