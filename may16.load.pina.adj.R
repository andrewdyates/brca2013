PINA.ADJ <- read.table("~/pina_adj_matrix.tab", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE, na.strings="")
save(PINA.ADJ, file="~/pina.adj.RData")
