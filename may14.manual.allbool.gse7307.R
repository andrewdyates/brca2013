library(Biobase)
load("../E7307.may14.genelvl.exprfilt.RData")
source("~/pymod/boolean_implication_fit/bool.R")
source("~/pymod/boolean_implication_fit/step.up.R")
source("~/pymod/boolean_implication_fit/all.pairs.R")

STEPS <- all.steps(M, do.plot=F)
CLS <- all.pairs.cls(M, STEPS, b)
save(STEPS, CLS, file="../GSE7307.may14.BOOL.RData")
