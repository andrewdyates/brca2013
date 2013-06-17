# THIS IS WRONG AND NEEDS TO USE FILTERED COPIES OF ALL BOOLEAN DISTANCES

python ~/code/qsub/script.py script="time python ~/code/all_bool_dist/script.py fname=~/brca/jun6.GSE31448.select.tab.b0.0968.z3.00.r0.20.err0.10.bool.tab" jobname="31448.TF.BOOL.D" n_ppn=8 email=True hours=6

python ~/code/qsub/script.py script="time python ~/code/all_bool_dist/script.py fname=~/brca/jun6.GSE7307.select.tab.b0.1541.z3.00.r0.10.err0.10.bool.tab" jobname="7307.TF.BOOL.D" n_ppn=8 email=True hours=6

