# coexpression

## Network density analysis to detect coexpression within genome-wide data

This script is the development version of the network density analysis 
used in our paper. In future, we hope to release a user-friendly, 
production version of the script, but for now we have shared the 
version of the code that we used for the analyses in the paper. 
It is far from perfect, but if you follow these instructions you 
should be able to replicate our findings, and generate your own. 
Please contact us through coexpression.net or GitHub if you have any 
problems.

Please refer to the supplementary methods for (http://biorxiv.org/content/early/2016/12/20/095349) for a full explanation of the method used here.

## TO DO

Add in additional info from https://raw.githubusercontent.com/baillielab/coexpression/master/README.md

## How to install

### Requirements

The coexpression code is mainly written in Python 2.7. Some code is written in C and OpenMP. In addition to Python, you will need:

* gcc (4.8.x recommended)
* [bedtools](https://bedtools.readthedocs.io) (2.x recommended)


### Python dependencies

An [Anaconda Python environment](https://www.anaconda.com/download) is recommended. The `environment.yml` file contains the Python dependencies. 

Create a conda enviroment:
```
$ conda env create -f environment.yml
$ source activate coexpression
```

## Example usage

First, run the prepare script:
```
$ python 0-prepare-coex.py -po -n 100 -w 2000 -p backcirc -f ${INPUT_FILE} -b ${BACKGROUND_FILE} 
```
this will create a directory called `../results-coex/BLAH`. Then for each of the permutations run:

```
$ RESULTS=../results-coex/BLAH
$ python 1-make-network.py -mf ${RESULTS}/permutation_store/f5ep300_100.NNNNNN.bed -wd ${RESULTS}
```

Finally, run the collate script:
```
$ python 2-collate-results.py -wd ${RESULTS}
```