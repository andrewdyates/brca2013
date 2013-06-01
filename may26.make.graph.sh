CLS="/Users/z/Dropbox/biostat/brca/may26.gse7307.k235.clusters/may26.k235.avglink.5pct.extant.collapsed.CLS.tab"
DCOR="/Users/z/Dropbox/biostat/brca/may26.gse7307.k235.clusters/may26.k235.avglink.5pct.extant.collapsed.DCOR.tab"
WEAK="/Users/z/Dropbox/biostat/brca/may26.gse7307.k235.clusters/may26.k235.avglink.5pct.extant.collapsed.WEAK.tab"

# BRCA1 in cluster 3
# BRCA2 in cluster 9

python /Users/z/Dropbox/biostat/git_repos/directed_graph/script.py min_d=0.40 cls_fname=$CLS dcor_fname=$DCOR graphviz_cmd=fdp outpath_prefix=~/Desktop/alltrans_brca_0.40_fdp


python /Users/z/Dropbox/biostat/git_repos/directed_graph/script.py min_d=0.5 cls_fname=$CLS dcor_fname=$DCOR graphviz_cmd=dot outpath_prefix=~/Desktop/alltrans_brca_0.5_dot_spline spline=true

python /Users/z/Dropbox/biostat/git_repos/directed_graph/script.py min_d=0.4 cls_fname=$CLS dcor_fname=$DCOR graphviz_cmd=dot outpath_prefix=~/Desktop/alltrans_brca_0.4_dot_spline spline=true


python /Users/z/Dropbox/biostat/git_repos/directed_graph/script.py min_d=0.5 cls_fname=$CLS dcor_fname=$DCOR graphviz_cmd=dot outpath_prefix=~/Desktop/alltrans_brca_0.5_dot_spline spline=true

python /Users/z/Dropbox/biostat/git_repos/directed_graph/script.py min_d=0.5 cls_fname=$CLS dcor_fname=$DCOR graphviz_cmd=dot outpath_prefix=~/Desktop/alltrans_brca_0.5_dot_ortho spline=ortho
