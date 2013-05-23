PINA.ADJ <- as.matrix(read.table("~/pina/pina_compiled_may22_2013_adjm.tab", sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE, na.strings="", check.names=F))
save(PINA.ADJ, file="~/pina/pina.may22.2013.adj.RData")
