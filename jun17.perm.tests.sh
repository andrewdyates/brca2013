python $HOME/code/qsub/script.py script="time python $HOME/code/permutation_test/script_single.py fname=$HOME/brca/jun17.R.GSE7307.TF.M.tab n=100000 dep=dcor > $HOME/brca/jun17.R.GSE7307.TF.M.tab.perm100000.out" n_ppn=8 hours=99 jobname=GSE7307_perm email=True

python $HOME/code/qsub/script.py script="time python $HOME/code/permutation_test/script_single.py fname=$HOME/brca/jun17.R.GSE31448.TF.M.tab n=100000 dep=dcor > $HOME/brca/jun17.R.GSE31448.TF.M.tab.perm100000.out" n_ppn=8 hours=99 jobname=GSE31448_perm email=True

python $HOME/code/qsub/script.py script="time python $HOME/code/permutation_test/script_single.py fname=$HOME/celegans/jun5.GSE2180.SCAN.select.tab.M.TFgold.tab n=100000 dep=dcor > $HOME/celegans/jun5.GSE2180.SCAN.select.tab.M.TFgold.tab.perm100000.out" n_ppn=8 hours=99 jobname=GSE2180_perm email=True


python $HOME/code/qsub/script.py script="time python $HOME/code/permutation_test/script_single.py fname=$HOME/brca/jun6.GSE7307.select.tab n=100000 dep=dcor > $HOME/brca/jun6.GSE7307.select.tab.perm100000.out" n_ppn=8 hours=99 jobname=GSE7307_all_perm email=True

python $HOME/code/qsub/script.py script="time python $HOME/code/permutation_test/script_single.py fname=$HOME/brca/jun6.GSE31448.select.tab n=100000 dep=dcor > $HOME/brca/jun6.GSE31448.select.tab.perm100000.out" n_ppn=8 hours=99 jobname=GSE31448_all_perm email=True

python $HOME/code/qsub/script.py script="time python $HOME/code/permutation_test/script_single.py fname=$HOME/celegans/jun5.GSE2180.SCAN.select.tab n=100000 dep=dcor > $HOME/celegans/jun5.GSE2180.SCAN.select.tab.perm100000.out" n_ppn=8 hours=99 jobname=GSE2180_all_perm email=True
