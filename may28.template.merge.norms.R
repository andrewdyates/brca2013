#"~/brca/GSE31448/normed"
library(Biobase)
load("~/brca/GSE31448/normed/SCAN.batch.1.RData")
load("SCAN.batch.2.RData")
load("SCAN.batch.3.RData")
load("UPC.batch.1.RData")
load("UPC.batch.2.RData")
load("UPC.batch.3.RData")

S.GSE31448 <- cbind(exprs(E.SCAN.batch.1), exprs(E.SCAN.batch.2), exprs(E.SCAN.batch.3))
U.GSE31448 <- cbind(exprs(E.UPC.batch.1), exprs(E.UPC.batch.2), exprs(E.UPC.batch.3))
save(S.GSE31448, U.GSE31448, file="~/brca/GSE31448/normed/GSE31448.SCANUPC.RData")
