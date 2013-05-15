library(Biobase)
load("../E7307.may14.genelvl.exprfilt.RData")
source("~/pymod/boolean_implication_fit/bool.R")
source("~/pymod/boolean_implication_fit/step.up.R")
source("~/pymod/boolean_implication_fit/all.pairs.R")
library(energy)

STEPS <- all.steps(M, do.plot=F)
save(STEPS, file="../STEPS.E7307.may14.R")

# compute select classes
# get BRCA1, BRCA2 indices
which(rownames(M) %in% c("BRCA1","BRCA2"))
[1]  2929 11098
BRCA1.CLS <- single.pairs.cls(M, STEPS, b, 2929)
BRCA2.CLS <- single.pairs.cls(M, STEPS, b, 11098)

summary(as.factor(BRCA1.CLS))
summary(as.factor(BRCA2.CLS))
# BRCA1
#   1    2    3    4    5 
# 610  142 2121 7375 3970
# BRCA2
#   1    2    3    4    5 
# 640   96 2007 7972 3503 

BRCA1.DCOR <- apply(M, 1, function(x) dcor(x,M[2929,]))
BRCA2.DCOR <- apply(M, 1, function(x) dcor(x,M[11098,]))
