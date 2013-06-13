python ~/code/qsub/script.py script="time python ~/code/all_bool_dist/script.py fname=~/brca/jun12.R.GSE31448.TF.BOOL.tab" options="#PBS -M yates.115.osu@gmail.com" jobname="31448.TF.BOOL" n_ppn=8 email=True walltime="6:00:00"
python ~/code/qsub/script.py script="time python ~/code/all_bool_dist/script.py fname=~/brca/jun12.R.GSE7307.TF.BOOL.tab" options="#PBS -M yates.115.osu@gmail.com" jobname="7307.TF.BOOL" n_ppn=8 email=True walltime="6:00:00"
python ~/code/qsub/script.py script="time python ~/code/all_bool_dist/script.py fname=~/brca/jun12.R.GSE31448.TF.WEAK.tab use_weak=True" options="#PBS -M yates.115.osu@gmail.com" jobname="31448.TF.WEAK" n_ppn=8 email=True walltime="6:00:00"
python ~/code/qsub/script.py script="time python ~/code/all_bool_dist/script.py fname=~/brca/jun12.R.GSE7307.TF.WEAK.tab use_weak=True" options="#PBS -M yates.115.osu@gmail.com" jobname="7307.TF.WEAK" n_ppn=8 email=True walltime="6:00:00"

