#!/bin/bash

# Input files and arguments
#
# Enter location of input and background .bed files; number of permutations and additional arguments
INPUT_FILE=/path/to/input.bed
BACKGROUND_FILE=/path/to/background.bed
N_PERMS=1
ARGUMENTS="-po -n ${N_PERMS} -p backcirc -f ${INPUT_FILE} -b ${BACKGROUND_FILE}"


# Activate conda environment
source activate coexpression

# Get location of results directory
RESULTS_DIR=$(python get_results_dir.py ${ARGUMENTS})
echo "Results directory is: ${RESULTS_DIR}"

# Prepare step
echo "Running prepare step"
python 0-prepare-coex.py ${ARGUMENTS} &> /dev/null

# Run all permutations
for i in $(seq 0 $N_PERMS); do
  echo "Running permutation: ${i}"
  python 1-make-network.py -mf ${RESULTS_DIR}/permutation_store/f5ep300_100.${i}.bed -wd ${RESULTS_DIR}
done

# Collate results
echo "Running collation step"
python 2-collate-results.py -wd ${RESULTS_DIR}
echo "Done"
