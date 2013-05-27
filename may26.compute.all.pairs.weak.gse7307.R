source("~/pymod/weak_implication_R/lib.R")
load("~/brca/E7307.may22.genelvl.exprfilt.RData")
GSE7307.WEAK <- all.pairs.weak(M)
save(GSE7307.WEAK, file="~/brca/E7307.may26.err4.th0.1.genelvl.exprfilt.WEAK.RData", err=4, th=0.1)
