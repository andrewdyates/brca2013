PINA.ADJ <- as.matrix(read.table("~/pina/pina_adj_matrix.tab", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE, na.strings="", check.names=F))
save(PINA.ADJ, file="~/pina/pina.adj.RData")
