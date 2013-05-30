source("~/pymod/weak_implication_R/lib.R")
load("~/brca/GSE31448.SCANUPC.RData")
GSE31448.WEAK <- all.pairs.weak(S.GSE31448, err=4, th=0.1)
save(GSE31448.WEAK, file="~/brca/E31448.may29.err4.th0.1.genelvl.exprfilt.WEAK.RData")
