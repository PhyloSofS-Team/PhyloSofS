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
[site](https://julialang.org/downloads/). Once, *Julia* is installed, open the [REPL](https://docs.julialang.org/en/v1/stdlib/REPL/) and install the needed packages by doing:

```julia
using Pkg
Pkg.add(["MIToS", "Statistics", "CSV", "Plots", "StatsPlots", "DataFrames", "BioAlignments", "BioStructures"])
```

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

To run the molecular modelling module you need the sequence database for
*HH-suite3* and the structural database (*PDB*) for *MODELLER*.
Follow the instructions in [docs/get_databases.md](shorturl.at/byIR8)
to set them up.

## Docker

To run [*PhyloSofS' Docker* image](https://cloud.docker.com/u/diegozea/repository/docker/diegozea/phylosofs)
you need to install *Docker* from the [*Docker* website](https://www.docker.com).

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

## Running PhyloSofS

You can run `phylosofs -h` to see the help and the list of arguments.

### 1. Phylogenetic Inference

```bash
phylosofs -P -s 100 --tree path_to_newick_tree --transcripts path_to_transcripts
```
### 2. Molecular modelling

```bash
phylosofs  -M -i path_to_input_files --hhlib path_to_hhsuite_folder --hhdb path_to_hblits_database(uniclust30)/basename --structdb path_to_hhpred(pdb70)/basename --allpdb path_to_the_cif_database --ncpu number_of_cpu --julia path_to_julia_executable
```
Please note that for the databases hhdb and structdb, you need to provide the path to the folder and also the base name of the files in it. For example, if the database uniclust30_2018_08 is located in /home, you need to write --hhdb /home/uniclust30_2018_08/uniclust30_2018_08

#### Docker example

```bash
phylosofs -M -i . --hhlib /app/hh-suite/ --hhdb /databases/uniclust30_2018_08/uniclust30_2018_08 --structdb /databases/pdb70_from_mmcif_latest/pdb70 --allpdb /databases/allpdb/ --ncpu 12
```

## Licence
The PhyloSofS package has been developed under the
[MIT License](https://github.com/PhyloSofS-Team/PhyloSofS/blob/master/LICENSE.txt).

## Contact
For questions, comments or suggestions feel free to contact
[Elodie Laine](mailto:elodie.laine@upmc.fr?subject=[GitHub]PhyloSofS) or
[Hugues Richard](mailto:hugues.richard@upmc.fr?subject=[GitHub]PhyloSofS)
