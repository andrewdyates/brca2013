library(Biobase)
load("../E7307.may14.genelvl.exprfilt.RData")
source("~/pymod/boolean_implication_fit/bool.R")
source("~/pymod/boolean_implication_fit/step.up.R")
source("~/pymod/boolean_implication_fit/all.pairs.R")
source("~/pymod/dependency_glyph_splom/lib.R")
library(energy)

STEPS <- all.steps(M, do.plot=F)
#save(STEPS, file="../STEPS.E7307.may14.R")


# compute select classes
# get BRCA1, BRCA2 indices
which(rownames(M) %in% c("BRCA1","BRCA2"))
#[1]  2929 11098
#2929 is BRCA2, 11098 is BRCA1
# i is row, y axis. 1 is YiX, 3 is XiY
BRCA2.CLS <- single.pairs.cls(M, STEPS, b, 2929)
BRCA1.CLS <- single.pairs.cls(M, STEPS, b, 11098)
save(BRCA1.CLS, BRCA2.CLS, file="../BRCA.CLS.RData")

summary(as.factor(BRCA1.CLS))
summary(as.factor(BRCA2.CLS))
# BRCA1
#   1    2    3    4    5 
# 610  142 2121 7375 3970
# BRCA2
#   1    2    3    4    5 
# 640   96 2007 7972 3503 

BRCA2.DCOR <- apply(M, 1, function(x) dcor(x,M[2929,]))
save(BRCA2.DCOR, file="../BRCA2.dcor.RData")
BRCA1.DCOR <- apply(M, 1, function(x) dcor(x,M[11098,]))
save(BRCA1.DCOR, file="../BRCA1.dcor.RData")

## TMP <- BRCA1.DCOR
## BRCA1.DCOR <- BRCA2.DCOR
## BRCA2.DCOR <- TMP

## TMP <- BRCA1.CLS
## BRCA1.CLS <- BRCA2.CLS
## BRCA2.CLS <- TMP
# ------
# May 15th checkpoint
# ------
load("../STEPS.E7307.may14.R")
load("../BRCA1.dcor.RData")
load("../BRCA2.dcor.RData")
load("../BRCA.CLS.RData")

pdf("../brca.expr.values.pdf")
plot(M[11098,],main="BRCA1")
plot(M[2929,],main="BRCA2")
dev.off()


pdf("../brca.dcor.hist.pdf")
hist(BRCA1.DCOR, main="BRCA1")
hist(BRCA2.DCOR, main="BRCA2")
dev.off()

# note: this splits gene lists separated by /// into names
qcm.gene.list <- as.character(read.table("../qcm_gene_list.txt")[,1])

D <- BRCA2.DCOR
CLS <- BRCA2.CLS
# top 100 genes by dcor
# top 30 genes by dcor per class
# boxplot per class
# class and dcor for select gene list
qq <- which(names(D) %in% qcm.gene.list)
length(qq)
#[1] 356
length(qcm.gene.list)
#[1] 426
 
 
brca.dcor.bool.plots <- function(D, CLS) {
  R <- list()
  R$topD <- sort(D, decreasing=TRUE)[1:100]
  R$clsh <- summary.plots.vector(CLS, D)
  qq <- which(names(D) %in% qcm.gene.list)
  R$qcm.dcor <- D[qq]
  R$qcm.cls <- CLS[qq]
  R
}

pdf("../summ.plots.brca1.pdf")
BRCA1.R <- brca.dcor.bool.plots(BRCA1.DCOR, BRCA1.CLS)
dev.off()

pdf("../summ.plots.brca2.pdf")
BRCA2.R <- brca.dcor.bool.plots(BRCA2.DCOR, BRCA2.CLS)
dev.off()

# table from BRCA1, BRCA2
BRCA.T <- data.frame(BRCA1.R$qcm.dcor, BRCA2.R$qcm.dcor, BRCA1.R$qcm.cls, BRCA2.R$qcm.cls)
write.table(BRCA.T, file="../gse7307.brcapush.targets.tab", sep="\t", quote=F, row.names=T, col.names=NA)

BRCA.TOP <- data.frame(names(BRCA1.R$topD), BRCA1.R$topD, names(BRCA2.R$topD), BRCA2.R$topD)
write.table(BRCA.TOP, file="../gse7307.brcapush.top.tab", sep="\t", quote=F, row.names=F, col.names=T)

# plots for targets only
pdf("../summ.plots.brca1.targets.pdf")
BRCA1.R <- brca.dcor.bool.plots(BRCA1.R$qcm.dcor, BRCA1.R$qcm.cls)
dev.off()

pdf("../summ.plots.brca2.targets.pdf")
BRCA2.R <- brca.dcor.bool.plots(BRCA2.R$qcm.dcor, BRCA2.R$qcm.cls)
dev.off()

## bool plots

qq <- which(rownames(M) %in% c(qcm.gene.list, "BRCA1", "BRCA2"))
M.Targets <- M[qq,]

pdf("../steps.plots.targets.pdf")
Targ.STEPS <- all.steps(M.Targets, do.plot=T)
dev.off()
pdf("../bool.plots.brca1.targets.pdf")
R <- single.pairs.cls(M.Targets, Targ.STEPS, b, which(rownames(M.Targets)=="BRCA1"), do.plot=TRUE)
dev.off()
pdf("../bool.plots.brca2.targets.pdf")
R <- single.pairs.cls(M.Targets, Targ.STEPS, b, which(rownames(M.Targets)=="BRCA2"), do.plot=TRUE)
dev.off()


