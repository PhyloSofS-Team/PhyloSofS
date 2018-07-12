# -*- coding: utf-8 -*-

# Copyright (c) 2014-2018: Adel Ait-hamlat, Elodie Laine, Lélia Polit
# and Diego Javier Zea.
# This code is part of the phylosofs package and governed by its license.
# Please see the LICENSE.txt file included as part of this package.

# module to reconstruct 3D models of isoforms

import glob
import os
import shutil
import subprocess
import warnings
import gzip
# import ConfigParser
# _config = ConfigParser.ConfigParser()

try:
    import modeller
    _modeller_message = ""
except ImportError:
    warnings.warn('modeller is not installed', ImportWarning)
    _modeller_message = "Please install modeller: " \
                        "https://salilab.org/modeller/"

# ------------------ UPLOAD INIT VALUES ----------------------------- #
#
#
# def init(configFile):
#
#    _config.read(os.path.expanduser(configFile))
#    # upload paths for librairies
#    HHBLITS = _config.get("PROGRAMS", "HHBLITS")
#    ADDSS = _config.get("PROGRAMS", "ADDSS")
#    HHMAKE = _config.get("PROGRAMS", "HHMAKE")
#    HHSEARCH = _config.get("PROGRAMS", "HHSEARCH")
#    HHMODEL = _config.get("PROGRAMS", "HHMODEL")
#    PROCHECK = _config.get("PROGRAMS", "PROCHECK")
#    NACCESS = _config.get("PROGRAMS", "NACCESS")
#    HHDB = _config.get("DATABASES", "HHDB")
#    STRUCTDB = _config.get("DATABASES", "STRUCTDB")
#    ALLPDB = _config.get("DATABASES", "ALLPDB")
#    NCPU = _config.get("OPTIONS", "NCPU")
#    CONTEXTLIB = _config.get("DATA", "CONTEXTLIB")
#
#    return(HHBLITS, ADDSS, HHMAKE, HHSEARCH, HHMODEL, PROCHECK, NACCESS, HHDB,
#           STRUCTDB, ALLPDB, NCPU, CONTEXTLIB)
#
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
        with open(seq[0] + ".fa", 'w') as f:
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
            s = 1
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
            seqs = readFastaMul(path, True)
            namedir = file.split('.')[0]
            if not os.path.exists(namedir):
                os.mkdir(namedir)
            os.chdir(namedir)
            writeFastas(seqs)
            ffiles = glob.glob('*.fa')
            for ff in ffiles:
                namertf = os.path.join(
                    transcriptsdir, ff.replace('.fa', '') + ".rtf")
                if os.path.isfile(namertf):
                    shutil.copy2(namertf, os.path.join(outputdir, namedir))
            os.chdir(outputdir)
    else:
        print transcriptsdir + " is not a folder (running in " + here + ")."
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
                print "can't read line"
            count += 1

        if count < nbLines:
            chain = ligne[21]
            ficOUT = code + "_" + chain + '.' + ext
            wfile2 = open(ficOUT, "w")
            wfile2.writelines(ligne)
            try:
                ligne = lines[count]
            except:  # TODO : Try not to use try here or use a better except.
                print "can't read line"
            count += 1
            while (ligne[21] == chain) and ((ligne.startswith("ATOM")) or
                                            (ligne.startswith("HETATM")) or
                                            (ligne.startswith("ANISOU"))) and \
                  (count < nbLines):
                wfile2.writelines(ligne)
                try:
                    ligne = lines[count]
                except:  # TODO : Try not to use try or use a better except.
                    print "can't read line"
                count += 1
            wfile2.close()


def assessNormalizedDopeScore(pdb):
    if _modeller_message != "":
        raise ImportError(_modeller_message)

    env = modeller.environ()
    env.libs.topology.read(file='$(LIB)/top_heav.lib')  # TODO: Test on Windows
    env.libs.parameters.read(file='$(LIB)/par.lib')
    # Read a model previously generated by Modeller's automodel class
    mdl = modeller.scripts.complete_pdb(env, pdb)
    zscore = mdl.assess_normalized_dope()
    return (zscore)


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


def model3D(fic, ALLPDB):
    if _modeller_message != "":
        raise ImportError(_modeller_message)

    modeller.log.verbose()  # request verbose output
    env = modeller.environ()  # create a new MODELLER environment to build ...
    # ... this model in

    seqs = readFastaMul(fic)
    seq = seqs[0][0].split('\n')[0]
    seq = seq.split(';')[1]

    knowns = []

    for i in range(1, len(seqs)):
        tmp = seqs[i][0].split('\n')[0]
        tmp = tmp.split(';')[1]
        knowns.append(tmp)
        kn = tmp.split('_')[0]
        base_name = 'pdb' + kn + '.ent'
        nam = base_name + '.gz'
        if not os.path.isfile('pdb' + kn + '.ent'):
            # shutil.copyfile(ALLPDB + nam, "./" + nam)
            shutil.copy2(ALLPDB + nam, nam)
            # os.system("gunzip ./" + nam)
            with gzip.open(nam, 'rb') as f_in:
                with open(base_name, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)

            splitChainsPDB('pdb' + kn + '.ent', kn, 'pdb')
    knowns = tuple(knowns)

    a = modeller.automodel.automodel(
        env, alnfile=fic, knowns=knowns, sequence=seq)
    a.starting_model = 1  # index of the first model
    a.ending_model = 1  # index of the last model
    # (determines how many models to calculate)
    cdir = os.getcwd()
    print("cdir = ", cdir)
    a.make()  # do the actual homology modeling

    return 1


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
    tmp = s.split("{\cf")
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
        myExons = readRTF(tmp + '.rtf')
        start = borders[2] - borders[0]
        end = len(myExons) - (borders[1] - borders[3] + 1)
        myExons = myExons[start:end]
        insertExons(myExons, tmp)
        res = 1
    except:
        print 'Error: could not annotate the 3D model for ' + tmp
        res = 0
    return res


def runModelProcess(HHBLITS, ADDSS, HHMAKE, HHSEARCH, HHMODEL, HHDB, STRUCTDB,
                    ALLPDB, NCPU, trans, selTemp, only3D, CONTEXTLIB):
    try:
        # tmp = trans[2:].split('.')[0]  # 'example.fa' --> 'ample'
        tmp = os.path.splitext(trans)[0]  # 'example.fa' --> 'example'
        if not only3D:
            # look for homologs and create the MSA alignment
            subprocess.call([
                HHBLITS, "-cpu", NCPU, "-i", trans, "-d", HHDB, "-oa3m",
                tmp + ".a3m", "-n", "1"
            ])
            # add secondary structure prediction to the MSA
            subprocess.call([ADDSS, tmp + ".a3m"])
            # generate a hidden Markov model (HMM) from the MSA
            subprocess.call([HHMAKE, "-i", tmp + ".a3m"])
            # search for homologs
            subprocess.call([
                HHSEARCH, "-cpu", NCPU, "-v", "1", "-i ", tmp + ".hhm", "-d",
                STRUCTDB, "-o", tmp + ".hhr", "-p", "20", "-P", "20", "-Z",
                "100", "-B", "100", "-seq", "1", "-aliw", "80", "-local",
                "-ssm", "2", "-norealign", "-sc", "1", "-dbstrlen", "10000",
                "-cs", CONTEXTLIB
            ])  # TODO : Is CONTEXTLIB defined?
        # create the alignment for modeller
        subprocess.call(
            [HHMODEL, "-i", tmp + ".hhr", "-pir", tmp + ".pir", "-m", selTemp])
        # treat the alignment file to remove N- and C-terminal loops
        borders = treatAli(tmp + '.pir')
        # generate the 3D models with Modeller
        model3D(tmp + '.pir', ALLPDB)
        annotate(trans, borders)
        res = 1
    except:  # TODO : Too general except
        print 'Error: could not build the 3D model for ' + tmp
        res = 0
    return res
