load("../jun17.BRCA.dep.matrices.TF.RData")
BR1e <- "672"; BR2e <- "675" # entrez gene IDs for select genes of interest

# filter any genes without at least a one-in-a-thousand relation
# ----------------------------------------
# GSE31448 perm test
# e-1:.12, e-2:.15, e-3:.18, e-4:.21, e-5:.25
D <- GSE31448.TF$DCOR
D[diag(dim(D)[1])==1] <- 0
GSE31448.tf.dcor.max <- apply(D, 1, max)
H <- hist(GSE31448.tf.dcor.max)
## $breaks
##  [1] 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90
## $counts
##  [1]   2  22  29  69 123 168 170 160 122  77  57  30  18   8
min(GSE31448.tf.dcor.max)
# [1] 0.211361
# Every gene is significantly related to something!


# ----------------------------------------
# GSE7307 perm test
# e-1:.08, e-2:.11, e-3:.13, e-4:.15, e-5:.18
D <- GSE7307.TF$DCOR
D[diag(dim(D)[1])==1] <- 0
GSE7307.tf.dcor.max <- apply(D, 1, max)
H <- hist(GSE7307.tf.dcor.max)
## $breaks
##  [1] 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 1.00
## $counts
##  [1]   1   5  16  27  58  95 135 144 124  91  65  33  18   2
min(GSE7307.tf.dcor.max)
# [1] 0.338931 # this is highly, highly significant

# map entrez IDs to symbols.
entrezSymTab <- read.table("~/code/ncbi_gene_info/human_gene_list.txt", header=F, sep="\t", as.is=T, quote="")
qq <- match(rownames(GSE31448.TF$DCOR), entrezSymTab[,1])
syms <- entrezSymTab[qq,2]
rownames(GSE31448.TF$DCOR) <- syms
colnames(GSE31448.TF$DCOR) <- syms
qq <- match(rownames(GSE7307.TF$DCOR), entrezSymTab[,1])
syms <- entrezSymTab[qq,2]
rownames(GSE7307.TF$DCOR) <- syms
colnames(GSE7307.TF$DCOR) <- syms

# replot using better fitted min and max
source("~/code/dependency_glyph_splom/lib.R")
GSE31448.BOOL.D.row <- as.dist(as.matrix(read.table("~/brca/jun17.R.GSE31448.TF.BOOL.tab.booldist.tab", header=T, row.names=1, check.names=F)))
stopifnot(all(rownames(GSE31448.BOOL.D.row)==colnames(GSE31448.BOOL.D.row)))
stopifnot(all(rownames(GSE31448.BOOL.D.row)==rownames(GSE31448.TF$BOOL)))
GSE7307.BOOL.D.row <- as.dist(as.matrix(read.table("~/brca/jun17.R.GSE7307.TF.BOOL.tab.booldist.tab", header=T, row.names=1, check.names=F)))
stopifnot(all(rownames(GSE7307.BOOL.D.row)==colnames(GSE7307.BOOL.D.row)))
stopifnot(all(rownames(GSE7307.BOOL.D.row)==rownames(GSE7307.TF$BOOL)))

# assign gene names to DCOR row


pdf("~/GSE7307.gsplom.tf.jun18.pdf", width=100, height=100)
R.GSE7307.TF <- splom(GSE7307.TF$BOOL, GSE7307.TF$DCOR, useRaster=T, sym=T, high.sat=T, MIN=0.13, MAX=0.8, row.cls.dist=GSE7307.BOOL.D.row)
dev.off()

pdf("~/GSE31448.gsplom.tf.jun18.pdf", width=100, height=100)
R.GSE31448.TF <- splom(GSE31448.TF$BOOL, GSE31448.TF$DCOR, useRaster=T, sym=T, high.sat=T, MIN=0.18, MAX=0.8, row.cls.dist=GSE31448.BOOL.D.row)
dev.off()

save(R.GSE7307.TF, R.GSE31448.TF, GSE31448.TF, GSE7307.TF, file="../jun18.BRCA.dep.matrices.TF.GSPLOM.RData")

# BRCA contradiction tests
# ----------------------------------------

GSE7307.sig <- GSE7307.TF$DCOR > 0.13
GSE31448.sig <- GSE31448.TF$DCOR > 0.18
sum(GSE7307.sig)
sum(GSE31448.sig)
sum(GSE31448.sig & !GSE7307.sig)
sum(!GSE31448.sig & GSE7307.sig)
