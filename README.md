# PhyloSofS

**A tool to model the evolution and structural impact of alternative splicing**

Status                     |Linux, OSX                 |Windows                    |Code Coverage
:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active) | [![Build Status](https://travis-ci.org/PhyloSofS-Team/PhyloSofS.svg?branch=master)](https://travis-ci.org/PhyloSofS-Team/PhyloSofS) | [![Build status](https://ci.appveyor.com/api/projects/status/jt1vvvawusokfx5c/branch/master?svg=true)](https://ci.appveyor.com/project/diegozea/phylosofs-fku85/branch/master) | [![codecov](https://codecov.io/gh/PhyloSofS-Team/PhyloSofS/branch/master/graph/badge.svg)](https://codecov.io/gh/PhyloSofS-Team/PhyloSofS) [![coveralls](https://coveralls.io/repos/github/PhyloSofS-Team/PhyloSofS/badge.svg?branch=master)](https://coveralls.io/github/PhyloSofS-Team/PhyloSofS?branch=master)

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
the **PhyloSofS' preprint** [available at *bioRxiv*](https://www.biorxiv.org/content/early/2017/03/23/119891).  

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

#### Optional external dependencies

- For the **phylogenetic inference**:  
  - [`PyGraphviz`](https://pygraphviz.github.io/) is needed to allow visualization of the phylogenies.
- For the **molecular modelling**:
  - [`HH-suite`](https://github.com/soedinglab/hh-suite) to identify homologous templates
  - [`MODELLER`](https://salilab.org/modeller/download_installation.html) to construct the 3D models.
  - [`NACCESS`](http://wolf.bms.umist.ac.uk/naccess/) to compute solvent accessible surface areas.
  - [`PROCHECK`](https://www.ebi.ac.uk/thornton-srv/software/PROCHECK/) to assess the quality of the models.

## Running PhyloSofS

### 1. Phylogenetic Inference

TO DO: Write usage examples

```
# ... to do ...
```

## Licence
The PhyloSofS package has been developed under the [MIT License](https://github.com/PhyloSofS-Team/PhyloSofS/blob/master/LICENSE.txt).

## Contact
For questions, comments or suggestions feel free to contact [Elodie Laine](mailto:elodie.laine@upmc.fr?subject=[GitHub]PhyloSofS) or [Hugues Richard](mailto:hugues.richard@upmc.fr?subject=[GitHub]PhyloSofS)
