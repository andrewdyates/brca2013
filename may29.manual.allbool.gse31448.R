library(Biobase)
load("../GSE31448.SCANUPC.RData")
source("~/pymod/boolean_implication_fit/bool.R")
source("~/pymod/boolean_implication_fit/step.up.R")
source("~/pymod/boolean_implication_fit/all.pairs.R")

M <- S.GSE31448
row.sd <- apply(M,1,sd)
b <- quantile(row.sd, 0.03)*2
## #0.09676083
## quantile(row.sd, 0.25)
## #0.1100014 (ok, looks reasonable)
## brca1.i <- which(rownames(M)=="672")
## brca2.i <- which(rownames(M)=="675")
## row.sd[brca1.i]
## # BRCA1 
## #      672 
## #0.2309285 
## #row.sd[brca2.i]
## #      675 
## #0.2158731
## sum(row.sd < row.sd[brca1.i]) / length(row.sd)
## #[1] 0.5912362

STEPS <- all.steps(M, do.plot=F)
CLS <- all.pairs.cls(M, STEPS, b)
save(STEPS, CLS, b, file="../GSE31448.may29.BOOL.RData")
