# coexpression

## Network density analysis to detect coexpression within genome-wide data

This script is the development version of the network density analysis 
used in our paper. In future, we hope to release a user-friendly, 
production version of the script, but for now we have shared the 
version of the code that we used for the analyses in the paper. 
It is far from perfect, but if you follow these instructions you 
should be able to replicate our findings, and generate your own. 
Please contact us through http://coexpression.net or GitHub if you have any 
problems.

Please refer to the supplementary methods for (http://biorxiv.org/content/early/2016/12/20/095349) for a full explanation of the method used here.

---
## Installation

### Requirements

The coexpression code is mainly written in Python 2.7. Some code is written in C and OpenMP. In addition to Python, you will need:

* gcc (4.8.x recommended)
* [bedtools](https://bedtools.readthedocs.io) (2.x recommended)
* supplementary files (see below)

### Compile C code

Compile the C code with `gcc`:
```
$ gcc -shared  -O3 -fPIC -fopenmp coexpression_v2.c -o coexpression_v2.so
```

### Python dependencies

An [Anaconda Python environment](https://www.anaconda.com/download) is recommended. The `environment.yml` file contains the Python dependencies. 

To create a conda enviroment, called 'coexpression':
```
$ conda env create -f environment.yml
$ source activate coexpression
```

### Configuration

The file `app.cfg` contains configuration information. This file needs to be changed to point at your copy of bedtools:

```
[directorypaths]
sourcefilesdir = ../supfiles-coex/
resdir = ../results-coex/
pathtobedtools = /path/to/bedtools
f5resource = http://fantom.gsc.riken.jp/zenbu/gLyphs/#config=ne92nJ20PhPv5ziW90qnND;loc=hg19::
```

### Supplementary files

Supplementary files should be downloaded from https://coexpression.roslin.ed.ac.uk/supportingfiles/supportingfiles.tar.gz 
and gunzipped into the directory `../supfiles-coex`. (The location of this directory can be changed in `app.cfg`)


### Test the installation

The following command runs a simple analysis and can be used to check that all the components and dependencies of 
coexpression are correctly installed.
```
$ python 0-prepare-coex.py -f ../supfiles-coex/test.bed -n 3 -x ../supfiles-coex/pure_f5p_condense_ptt_minimal_av.expression -v
```
Output files will be written to the directory: `../results-coex/test_complete_CIRCULAR_pj0.1_pure_f5p_condense_ptt_minimal_av`

---
## Options

The prepare script, `0-prepare-coex.py` takes the following options:

* `-f` [REQUIRED] specify the input file (modified bed format) containing the test set of interest - that is, the genomic locations of SNPs associated with a given phenotype
* `-p` the permutation mode: 'circular', 'post', 'backcirc' or 'background'. Default is 'circular'
* `-b` path to background file (modified bed format) detailing ALL variants genotypes in a particular study 
* `-n` number of permutations to run - default is 100
* `-t` number of threads to use 
* `-x` specify an expression data file (table of expression values for each feature - in the default mode, features are FANTOM5 TSS) - must match up to feature bed file
* `-q` specify an feature bed data file (table of genomic locations for each feature - in the default mode, features are FANTOM5 TSS) - must match up to expression data file
* `-w` a window size in base pairs around each TSS (or other input region). SNPs mapping within this window are taken to label a given TSS (or other input region) as being putatively phenotye-associated.
* `-j` specify the empirical p-value for correlation at which nearby promoters/enhancers are taken to belong to the same group (see supplementary methods) 
* `-e` user@example.com user email address for notification 
* `-v` verbose mode 
* `-a` include anticorrelations 
* `-i` don't use iterative removal of the most significantly-coexpressed regions (carries risk of false-positives) 
* `-z` don't use existing permutations already calculated on previous runs of this script
* `-po` prepare only: stop after preparation step to allow each permutation to be run separately (and possibly in parallel)
* `-cm` correlation measure: 'Spearman' or 'Pearson'. Default is 'Spearman'
* `-s` precision of node-specific p-value. Default is 1 (perfect), values closer to zero will save time. 

The scripts `1-make-network.py` and `2-collate-results.py` take options:

* `-wd` working directory: specify the path to the output directory created by `0-prepare-coex.py`. 
Options given to `0-prepare-coex.py` are read in from the directory 
* `-mf` path to map-file

---
## Example usage

First, run the prepare script:
```
$ python 0-prepare-coex.py -po -n 100 -p backcirc -f <inputfile.bed> -b <backgroundfile.bed> 
```
This will create an output directory in `../results-coex/` (configurable in `app.cfg`), named according to the input 
filenames and command-line arguments used. 

The prepare script creates one `.bed` file for each permutation (the number of permutations is specified by the `-n` 
option).

Each permutation can be run independently, e.g. to run the first (n=0):
```
$ RESULTS=../results-coex/name_of_results_directory
$ N=0
$ python 1-make-network.py -mf ${RESULTS}/permutation_store/f5ep300_100.$N.bed -wd ${RESULTS}
```

Finally, run the collate script:
```
$ python 2-collate-results.py -wd ${RESULTS}
```