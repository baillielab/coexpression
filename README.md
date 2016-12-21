# coexpression

===Network density analysis to detect coexpression within genome-wide data ===

This script is the development version of the network density analysis used in our paper. In future, we hope to release a user-friendly, production version of the script, but for now we have shared the version of the code that we used for the analyses in the paper. It is far from perfect, but if you follow these instructions you should be able to replicate our findings, and generate your own. Please contact us through coexpression.net or github if you have any problems.

Please refer to the supplementary methods for http://biorxiv.org/content/early/2016/12/20/095349 for a full explanation of the method used here. 

==1. Installation ==
Tested on scientific linux 6 using python 2.7
Dependent on the following python packages:
- getopt
- sys
- os
- time
- threading
- resource
- gc
- datetime
- subprocess 
- Popen
- PIPE
- STDOUT
- ctypes 
- string
- copy
- stat
- random
- math
- inspect
- itertools
- operator
- bisect 
- scipy 
- pandas
- networkx
- timeit
- json
- numpy

Copy the "run_coexpression_vX.py" script into an empty directory. Download the supporting files here:
https://coexpression.roslin.ed.ac.uk/supportingfiles/supportingfiles.tar.gz
and place them in a subdirectory called "supporting_files")


==2. Running the script ==

Example:

python run_coexpression_v1.0.py -f supportingfiles/test.bed -n 100

Options:
-f <filepath> specify the input file (modified bed format) containing the test set of interest
-b <filepath> path to background file (modified bed format) detailing ALL variants genotypes in a particular study
-n <integer> number of permutations to run
-t <integer> number of threads to use
-x <filepath> specify the expression data file (table of expression values for each FANTOM5 promoter/enhancer)
-j <float> specify the empirical p-value for correlation at which nearby promoters/enhancers are taken to belong to the same group (see supplementary methods)
-e <user@example.com> user email address for notification
-g use circular permutations on the background file (recommended)
-p use post-mapping permutations only (less powerful, eliminates enrichment signal)
-v verbose mode
-a include anticorrelations
-i don't use iterative removal of the most significantly-coexpressed regions (carries risk of false-positives)
-z don't use existing permutations already calculated on previous runs of this script. 


==Modified bed format==

Tab-seperated file with the following columns:
chr	start_position	end_position	[optional]snp_id	[optional]p_value
