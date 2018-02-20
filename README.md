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
## TO DO

Add in additional info from https://raw.githubusercontent.com/baillielab/coexpression/master/README.md

---
## How to install

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
and gunzipped into the directory `../supfiles-coex`. (The location of this directory can be changed in `app.cfg`.)

---

## Example usage

First, run the prepare script:
```
$ python 0-prepare-coex.py -po -n 100 -w 300 -p backcirc -f ${INPUT_FILE} -b ${BACKGROUND_FILE} 
```
this will create an output directory in `../results-coex/` (configurable in `app.cfg`), named according to the input filenames and command-line arguments 
used. 

The prepare script creates one `.bed` file for each permutation (the number of permutations is specified by the `-n` option).

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