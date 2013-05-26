library(Biobase)
load("~/brca/E7307.may14.genelvl.exprfilt.nice.pkl.DCOR.tab.RData")
DCOR <- M
load("~/brca/GSE7307.may14.BOOL.RData")
load("~/tftargets/tf.adj.may22.2013.RData")
load("~/pina/pina.may22.2013.adj.RData")
load("../E7307.may22.genelvl.exprfilt.RData")
exprM <- M
source("~/pymod/dependency_glyph_splom/lib.R")
qcm <- readLines("~/pymod/brca_qcm_list/qcm_gene_list_clean.txt")
select.syms <- readLines("BRCA.network.gene.list.txt") # 1891
sym.remap <- read.table("E7307.genesym.remap.txt", stringsAsFactors=F) # [1] 14218     2

# First pass summary plots
pdf("gsplom.summary.all.pdf")
summary.plots(CLS, DCOR, sym=T)
dev.off()

# Load weaks
load("../E7307.may22.genelvl.exprfilt.WEAK.RData")
# Load PCC
load("../E7307.may14.genelvl.exprfilt.nice.pkl.PEARSON.values.tab.RData")
PCC <- M

# CLS, DCOR, PCC, GSE7307.WEAK, STEPS, E, exprM

# Verify ID mappings
all(rownames(CLS)==colnames(CLS)) # T
all(rownames(CLS)==rownames(DCOR)) # T
all(rownames(CLS)==rownames(GSE7307.WEAK)) # FALSE
all(rownames(exprM)==rownames(GSE7307.WEAK)) # T
all(colnames(GSE7307.WEAK)==rownames(GSE7307.WEAK)) # T
all(rownames(exprM)==colnames(GSE7307.WEAK)) # T
all(names(STEPS)==rownames(exprM)) # FALSE
all(names(STEPS)==rownames(CLS)) # T
# OK: exprM and WEAK are updated, DCOR, CLS, and STEPS are not.
rownames(DCOR) <- rownames(exprM)
colnames(DCOR) <- rownames(exprM)

rownames(CLS) <- rownames(exprM)
colnames(CLS) <- rownames(exprM)

rownames(PCC) <- rownames(exprM)
colnames(PCC) <- rownames(exprM)

names(STEPS) <- rownames(exprM)
# Package up everything into a single file for easy loading later.
WEAK <- GSE7307.WEAK
save(CLS, DCOR, WEAK, exprM, E, STEPS, b, PCC, sym.remap, select.syms, file="../GSE7307.data.may25.RData")

# select select syms from each matrix
qq <- rownames(CLS) %in% select.syms
length(select.syms) # 1891
sum(qq) # 1347

# extract them from each matrix
CLS.s <- CLS[qq,qq]
PCC.s <- PCC[qq,qq]
DCOR.s <- DCOR[qq,qq]
WEAK.s <- WEAK[qq,qq]
exprM.s <- exprM[qq,]
STEPS.s <- STEPS[qq]
# save a file of only that selection
save(CLS.s,DCOR.s,WEAK.s,PCC.s,exprM.s,STEPS.s,b, file="../GSE7307.select.may25.RData")
writeLines(rownames(CLS.s), "../selected_symbols_w_data.txt")
