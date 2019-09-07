# PhyloSofS

**A tool to model the evolution and structural impact of alternative splicing**

Status                     |Linux, OSX                 |Windows                    
:-------------------------:|:-------------------------:|:-------------------------
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) | [![Build Status](https://travis-ci.org/PhyloSofS-Team/PhyloSofS.svg?branch=master)](https://travis-ci.org/PhyloSofS-Team/PhyloSofS) | [![Build status](https://ci.appveyor.com/api/projects/status/jt1vvvawusokfx5c/branch/master?svg=true)](https://ci.appveyor.com/project/diegozea/phylosofs-fku85/branch/master)

*PhyloSofS* (**Phylo**genies of **S**plicing is**of**orms **S**tructures) is a
fully automated computational tool that infers plausible evolutionary scenarios
explaining a set of transcripts observed in several species and models the
three-dimensional structures of the produced protein isoforms.  
The phylogenetic reconstruction algorithm relies on a combinatorial approach
and the maximum parsimony principle. The generation of the isoforms' 3D models
is performed using comparative modeling.  

#### Case study

PhyloSofS was applied to the c-Jun N-terminal kinase (JNK) family
(60 transcripts in 7 species). It enabled to date the appearance of an
alternative splicing event (ASE) resulting in substrate affinity modulation in
the ancestor common to mammals, amphibians and fishes, and to identify key
residues responsible for such modulation. It also highlighted a new ASE
inducing a large deletion, yet conserved across several species. The resulting
isoform is stable in solution and could play a role in the cell. More details
about this case study, together with the algorithm description, can be found in
the **PhyloSofS' preprint**
[available at *bioRxiv*](https://doi.org/10.1101/119891).  

## Installation

### 1. Download

You can clone this PhyloSofS package using [`git`](https://git-scm.com/):

```
git clone https://github.com/PhyloSofS-Team/PhyloSofS.git
```

### 2. Install

Then, you can access the cloned `PhyloSofS` folder and install the package
using Python 3's [`pip`](https://pip.pypa.io/en/stable/installing/):

```
cd PhyloSofS
python -m pip install .
```

### 3. Install dependencies

#### Phylogenetic inference

To run the phylogenetic module of *PhyloSofS*, you need to have
[*Graphviz*](https://graphviz.org/) installed.  
The easiest way to install *Graphviz* in...
  - **Debian/Ubuntu** is: `sudo apt-get install graphviz`
  - **Windows** is using [*Chocolatey*](https://chocolatey.org/): `choco install graphviz`
  - **macOS** is using [*Homebrew*](https://brew.sh/index): `brew install graphviz`

#### Molecular modelling

The molecular modelling pipeline depends on *Julia*, *HH-suite3* and *MODELLER*.
This module can only run on *Unix* systems (because of the *HH-suite*). To
alleviate that, we offer a *Docker* image with all these dependencies installed
(see the *Docker* section for more details).

##### Julia

You can download *Julia* 1.1.1 binaries from its
[site](https://julialang.org/downloads/).

###### LibZ

Some *BioJulia* packages can need *LibZ* to precompile. If you found a related
error, you can install *LibZ* from its [site](http://zlib.net/).
In *Ubuntu 18.04* you can install it by doing: `sudo apt-get install zlib1g-dev`

##### HH-suite3

Clone our *HH-suite* fork at
[AntoineLabeeuw/hh-suite](https://github.com/AntoineLabeeuw/hh-suite)
and follow the
[*Compilation*](https://github.com/AntoineLabeeuw/hh-suite#compilation)
instructions in its *README.md* file.

##### MODELLER

PhyloSofS needs MODELLER version 9.21. Follow the instructions in the
[*MODELLER* site](https://salilab.org/modeller/download_installation.html)
to install it and get the license key.

##### Databases

To run the molecular modelling module you need the
[*HH-suite* databases](http://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/):
  - Sequence database: `uniclust30_yyyy_mm_hhsuite.tar.gz`
  (we have tested PhyloSofS using `20180_08` as `yyyy_mm`)
  - Structural database: `pdb70_from_mmcif_latest.tar.gz`

The needed *mmCIF* *PDB* files for *MODELLER* are downloaded on demand, if
there are not present, in an indicated folder.

To set up the databases, you can use the script `setup_databases` (recommended).
Alternatively, a manual installation can be performed following the instructions
in [docs/get_databases.md](https://github.com/PhyloSofS-Team/PhyloSofS/blob/master/doc/get_databases.md).

###### Using the `setup_databases` script

The `setup_databases` downloads and decompress the needed databases. It creates
the following folder structure that can be easily used by *PhyloSofS* with the
`--databases` argument:

```
databases
 ├── pdb
 ├── pdb70
 └── uniclust
```

You can do `setup_databases -h` to know more about the script and its arguments.

## Docker

To run [*PhyloSofS' Docker* image](https://cloud.docker.com/u/diegozea/repository/docker/diegozea/phylosofs)
you need to install [*Docker*](https://www.docker.com) following [these instructions](https://docs.docker.com/v17.09/engine/installation/).

The following example is going to run *PhyloSofS' Docker* image using
*Windows PowerShell*. Databases for the molecular modelling module stored in
`D:\databases` are going to be mounted in `/databases` and the local directory
in `/project`. The actual folder is `${PWD}` in *Windows PowerShell*, `%cd%` in
*Windows Command Line* (*cmd*), and `$(pwd)` in *Unix*.

```
docker run -ti --rm --mount type=bind,source=d:\databases,target=/databases --mount type=bind,source=${PWD},target=/project diegozea/phylosofs
```

After this, we have access to the `bash` terminal of an *Ubuntu 18.04* image
with *PhyloSofS* and all its dependencies installed. You only need to indicate
your [MODELLER license key](https://salilab.org/modeller/registration.html) to
use *PhyloSofS*. To do that, you run the following command after replacing
`license_key` with your *MODELLER* license key:

```bash
sed -i 's/xxx/license_key/' /usr/lib/modeller9.21/modlib/modeller/config.py
```

### Installing databases with *Windows* as a host

If you are using the *PhyloSofS' Docker* image, you must know that errors can
occur when very large files are being written to bind-mounted *NTFS* file
systems. This happens particularly when `setup_databases` is run because it
tries to download and decompress large files. To avoid this problem, you can
install *PhyloSofS* on *Windows* and run `setup_databases.exe` to set up the
databases before using the docker image.  

## Running PhyloSofS

You can run `phylosofs -h` to see the help and the list of arguments.

### 1. Phylogenetic Inference

```bash
phylosofs -P -s 100 --tree path_to_newick_tree --transcripts path_to_transcripts
```
### 2. Molecular modelling

If databases where installed using `setup_databases` and the *HH-Suite3*
scripts and programs are in the executable paths, then you can run:

```bash
phylosofs -M -i path_to_input_dir --databases path_to_databases_folder
```

*PhyloSofS* is going to look for `transcripts.pir` files in the folder and
sub-folders of `path_to_input_dir` to perform the homology modelling of each
sequence in those files.

If you have a more manual installation of the databases and/or the *HH-Suite3*
scripts and programs are not in the path:

```bash
phylosofs -M -i path_to_input_dir --hhlib path_to_hhsuite_folder --hhdb path_to_uniclust_database/uniclust_basename --structdb path_to_pdb70/pdb70 --allpdb path_mmcif_pdb_cache_folder
```

Please note that for the databases `--hhdb` and `--structdb`, you need to
provide the path to the folder and also the basename of the files in it.
For example, if the database `uniclust30_2018_08` is located in `/home`, you
need to write:  

```
--hhdb /home/uniclust30_2018_08/uniclust30_2018_08
```

You can also find useful the arguments:
 - `--ncpu number_of_cpu`
 - `--julia path_to_julia_executable`

#### Docker example

If you installed the databases using `setup_databases` in a folder that you
have bind-mounted to `/databases`, then you only need to run:

```bash
phylosofs -M --databases /databases
```

Because the *Docker* image has *HH-Suite3* installed with its programs and
scripts in the executable paths.

## Licence
The PhyloSofS package has been developed under the
[MIT License](https://github.com/PhyloSofS-Team/PhyloSofS/blob/master/LICENSE.txt).

## Contact
For questions, comments or suggestions feel free to contact
[Elodie Laine](mailto:elodie.laine@upmc.fr?subject=[GitHub]PhyloSofS) or
[Hugues Richard](mailto:hugues.richard@upmc.fr?subject=[GitHub]PhyloSofS)
