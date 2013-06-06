GSE31448=$HOME/brca/june5_data/GSE31448.SCAN.tab
GSE7307=$HOME/brca/june5_data/GSE7307.SCAN.tab

python $HOME/pymod/dependency_matrix/compile_to_R_dispatch_wrapper.py fname=$GSE31448 computers=[\"PCC\",\"Dcor\"] outdir=$HOME/brca/june5_data/GSE31448_dep n_nodes=10 n_ppn=12 hours=72

python $HOME/pymod/dependency_matrix/compile_to_R_dispatch_wrapper.py fname=$GSE7307 computers=[\"PCC\",\"Dcor\"] outdir=$HOME/brca/june5_data/GSE7307_dep n_nodes=10 n_ppn=12 hours=72

#python $HOME/pymod/boolean_implication_fit_py/script_qsub.py fname=$GSE31448 d=bool
#python $HOME/pymod/boolean_implication_fit_py/script_qsub.py fname=$GSE31448 d=weak
#python $HOME/pymod/boolean_implication_fit_py/script_qsub.py fname=$GSE7307 d=bool
#python $HOME/pymod/boolean_implication_fit_py/script_qsub.py fname=$GSE7307 d=weak
