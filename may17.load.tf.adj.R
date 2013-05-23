TF.ADJ <- as.matrix(read.table("~/tftargets/all_tf_with_targs_may_21_2013_adjm.tab", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE, na.strings="", check.names=F))
save(TF.ADJ, file="~/tftargets/tf.adj.may22.2013.RData")
