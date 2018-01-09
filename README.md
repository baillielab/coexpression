==Installation==

=Dependencies=
The following python modules must be installed:
- networkx



==Example commands==
=1 Run a small test analysis=
- small input bed file
- run 3 permutations
- use a tiny expression file (a tenth of the size of the standard one)
python 0-run-coex.py -f ../supfiles-coex/test.bed -n 3 -x ../supfiles-coex/pure_f5p_condense_ptt_minimal_av.expression -v


=2 Run a single network analysis=
python 1-make-network.py -ef ../supfiles-coex/pure_f5p_condense_ptt_minimal_av.expression -fc /mnt/ris-fas1a/linux_groups2/baillie_grp/coexpression/coex-app/../supfiles-coex/f5ep300_100.bed -cf /mnt/ris-fas1a/linux_groups2/baillie_grp/coexpression/coex-app/../supfiles-coex/pure_f5p_condense_ptt_minimal_av.expression.spearmanlist -mf /mnt/ris-fas1a/linux_groups2/baillie_grp/coexpression/coex-app/../results-coex/test_complete_CIRCULAR_pj0.1_pure_f5p_condense_ptt_minimal_av/premapping_permutation_store/f5ep300_100_test.snps_mapped.0.bed -cm Spearman -sl Spearman_0_3CIRCULARperms -wd /mnt/ris-fas1a/linux_groups2/baillie_grp/coexpression/results-coex/test_complete_CIRCULAR_pj0.1_pure_f5p_condense_ptt_minimal_av -v 

python 1-make-network.py -ef ../supfiles-coex/pure_f5p_condense_ptt_minimal_av.expression -fc /mnt/ris-fas1a/linux_groups2/baillie_grp/coexpression/supfiles-coex/f5ep300_100.bed -cf /mnt/ris-fas1a/linux_groups2/baillie_grp/coexpression/coex-app/../supfiles-coex/pure_f5p_condense_ptt_minimal_av.expression.spearmanlist -mf /mnt/ris-fas1a/linux_groups2/baillie_grp/coexpression/coex-app/../results-coex/test_complete_CIRCULAR_pj0.1_pure_f5p_condense_ptt_minimal_av/permutation_store/f5ep300_100.0.bed -wd /mnt/ris-fas1a/linux_groups2/baillie_grp/coexpression/coex-app/../results-coex/test_complete_CIRCULAR_pj0.1_pure_f5p_condense_ptt_minimal_av -sl coex