TF.ADJ <- read.table("~/tftargets/tf_adj_matrix.tab", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE, na.strings="", check.names=F)
save(TF.ADJ, file="~/tftargets/tf.adj.RData")

TF.ADJ <- read.table("~/tftargets/tf_adj_matrix_aliasmapped.tab", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE, na.strings="", check.names=F)
save(TF.ADJ, file="~/tftargets/tf.adj.aliasmapped.RData")
