tf.entrez <- readLines("~/code/transcription_factors/all_tf_with_targs_may_21_2013.csv.entrez.tfonly.txt")
# load newly-generated entrez-gene level expression matrices
load("~/brca/jun6.GSE7307.GSE31448.select.RData")
## #"b.31448"    "b.7307"     "S.GSE31448" "S.GSE7307"
dim(S.GSE7307)
## #[1] 13766   677
dim(S.GSE31448)
## #[1] 14794   357

load("~/brca/jun6.GSE31448.select.tab.b0.0968.z3.00.r0.80.err0.10.bool.tab.RData")
GSE31448.BOOL <- M
load("~/brca/jun6.GSE31448.select.tab.err4.th0.2000.weak.tab.RData")
GSE31448.WEAK <- M
load("~/brca/GSE31448.SCAN.pkl.DCOR.values.tab.RData")
GSE31448.DCOR <- M

load("~/brca/jun6.GSE7307.select.tab.b0.1541.z3.00.r0.90.err0.10.bool.tab.RData")
GSE7307.BOOL <- M
load("~/brca/jun6.GSE7307.select.tab.err3.th0.2000.weak.tab.RData")
GSE7307.WEAK <- M
load("~/brca/GSE7307.SCAN.pkl.DCOR.values.tab.RData")
GSE7307.DCOR <- M

# select TF only
BR1e <- "672"; BR2e <- "675" # entrez gene IDs for select genes of interest
# genes per class for BR1, BR2
# all genes

select.and.align <- function(BOOL, DCOR, M, WEAK=NULL, idlist) {
    R = list()
    idlist <- c(idlist, BR1e, BR2e)
    qq <- rownames(BOOL) %in% idlist
    R$BOOL <- BOOL[qq,qq]
    qq <- rownames(DCOR) %in% idlist
    R$DCOR <- DCOR[qq,qq]
    qq <- rownames(M) %in% idlist
    R$M <- M[qq,]
    qq <- match(rownames(R$DCOR), rownames(R$BOOL))
    R$BOOL <- R$BOOL[na.omit(qq),na.omit(qq)]
    R$DCOR <- R$DCOR[!is.na(qq),!is.na(qq)]
    qq <- match(rownames(R$DCOR), rownames(R$M))
    stopifnot(all(!is.na(qq)))
    R$M[qq,]
    if(!is.null(WEAK)) {
        qq <- rownames(WEAK) %in% idlist
        R$WEAK <- WEAK[qq,qq]
        qq <- match(rownames(R$DCOR), rownames(R$WEAK))
        R$WEAK <- R$WEAK[na.omit(qq),na.omit(qq)]
        R$DCOR <- R$DCOR[!is.na(qq),!is.na(qq)]
        R$BOOL <- R$BOOL[!is.na(qq),!is.na(qq)]
    }
    stopifnot(all(rownames(R$DCOR)==colnames(R$DCOR)))
    stopifnot(all(rownames(R$BOOL)==colnames(R$BOOL)))
    stopifnot(all(rownames(R$BOOL)==rownames(R$DCOR)))
    if(!is.null(WEAK)) {
        stopifnot(all(rownames(R$WEAK)==colnames(R$WEAK)))
        stopifnot(all(rownames(R$WEAK)==rownames(R$DCOR)))
    }
    R
}

GSE31448.TF <- select.and.align(GSE31448.BOOL, GSE31448.DCOR, S.GSE31448, GSE31448.WEAK, tf.entrez)
GSE7307.TF <- select.and.align(GSE7307.BOOL, GSE7307.DCOR, S.GSE7307, GSE7307.WEAK, tf.entrez)
dim(GSE31448.TF$BOOL)
#[1] 1055 1055
dim(GSE7307.TF$BOOL)
#[1] 814 814
dim(S.GSE7307)
#[1] 13766   677
dim(S.GSE31448)
#[1] 14794   357
save(GSE31448.TF, GSE7307.TF, file="../jun17.BRCA.dep.matrices.TF.RData")
stopifnot(FALSE)

# EXAMINE BRCA* distributions
gse31448.br1i <- which(rownames(GSE31448.TF$BOOL)==BR1e)
gse31448.br2i <- which(rownames(GSE31448.TF$BOOL)==BR2e)
gse7307.br1i <- which(rownames(GSE7307.TF$BOOL)==BR1e)
gse7307.br2i <- which(rownames(GSE7307.TF$BOOL)==BR2e)
source("~/code/dependency_glyph_splom/lib.R")
pdf("~/GSE31448.brca1.summary.plots.jun17.pdf")
H <- summary.plots.vector(GSE31448.TF$BOOL[gse31448.br1i,], GSE31448.TF$DCOR[gse31448.br1i,])
dev.off()
pdf("~/GSE31448.brca2.summary.plots.jun17.pdf")
H <- summary.plots.vector(GSE31448.TF$BOOL[gse31448.br2i,], GSE31448.TF$DCOR[gse31448.br2i,])
dev.off()
pdf("~/GSE7307.brca1.summary.plots.jun17.pdf")
H <- summary.plots.vector(GSE7307.TF$BOOL[gse7307.br1i,], GSE7307.TF$DCOR[gse7307.br1i,])
dev.off()
pdf("~/GSE7307.brca2.summary.plots.jun17.pdf")
H <- summary.plots.vector(GSE7307.TF$BOOL[gse7307.br2i,], GSE7307.TF$DCOR[gse7307.br2i,])
dev.off()

b1.7 <- sort(GSE7307.TF$DCOR[gse7307.br1i,], decreasing=T)
b1.3 <- sort(GSE31448.TF$DCOR[gse31448.br1i,], decreasing=T)

length(names(b1.7[b1.7 > 0.3])) # 185
length(names(b1.3[b1.3 > 0.3])) # 26
shared.tf.b1 <- intersect(names(b1.7[b1.7 > 0.3]), names(b1.3[b1.3 > 0.3]))
length(shared.tf.b1) # 16

sort(b1.3[match(shared.tf.b1, names(b1.3))], decreasing=T)
sort(b1.7[match(shared.tf.b1, names(b1.7))], decreasing=T)

b1.7b <- GSE7307.TF$BOOL[gse7307.br1i,]
b1.3b <- GSE31448.TF$BOOL[gse31448.br1i,]
summary(as.factor(b1.7b))
#  0   1   2   3   4   5 
#125  26   9  89 379 186
summary(as.factor(b1.3b))
#  0   1   2   4   5 
#133   5   1 902  14 

b1.7b[match(shared.tf.b1, names(b1.7b))]
##    672   2305   2146 144455   7468   6941   1869  79733   4605   4603  25888 
##      2      2      2      1      2      1      1      3      2      1      1 
## 221504   3607   6925  26137   6945 
##      2      0      5      4      3 
b1.3b[match(shared.tf.b1, names(b1.3b))]
##    672   2305   2146 144455   7468   6941   1869  79733   4605   4603  25888 
##      2      4      4      4      4      4      4      4      4      4      4 
## 221504   3607   6925  26137   6945 
##      4      4      4      4      4 

b2.7 <- sort(GSE7307.TF$DCOR[gse7307.br2i,], decreasing=T)
b2.3 <- sort(GSE31448.TF$DCOR[gse31448.br2i,], decreasing=T)
length(names(b2.7[b2.7 > 0.3])) # 322
length(names(b2.3[b2.3 > 0.3])) # 176
shared.tf.b2 <- intersect(names(b2.7[b2.7 > 0.3]), names(b2.3[b2.3 > 0.3]))
length(shared.tf.b2) # 87

sort(b2.7[match(shared.tf.b2, names(b2.7))], decreasing=T)
sort(b2.3[match(shared.tf.b2, names(b2.3))], decreasing=T)

b2.7b <- GSE7307.TF$BOOL[gse7307.br2i,]
b2.3b <- GSE31448.TF$BOOL[gse31448.br2i,]
summary(as.factor(b2.7b))
summary(as.factor(b2.3b))
b2.7b[match(shared.tf.b2, names(b2.7b))]
b2.3b[match(shared.tf.b2, names(b2.3b))]

# ==================================================
# Transcription Factor Network
# ==================================================
# write matrices to file
write.table(GSE31448.TF$BOOL, file="../jun17.R.GSE31448.TF.BOOL.tab", sep="\t")
write.table(GSE31448.TF$WEAK, file="../jun17.R.GSE31448.TF.WEAK.tab", sep="\t")
write.table(GSE31448.TF$DCOR, file="../jun17.R.GSE31448.TF.DCOR.tab", sep="\t")
write.table(GSE31448.TF$M, file="../jun17.R.GSE31448.TF.M.tab", sep="\t")

write.table(GSE7307.TF$BOOL, file="../jun17.R.GSE7307.TF.BOOL.tab", sep="\t")
write.table(GSE7307.TF$WEAK, file="../jun17.R.GSE7307.TF.WEAK.tab", sep="\t")
write.table(GSE7307.TF$DCOR, file="../jun17.R.GSE7307.TF.DCOR.tab", sep="\t")
write.table(GSE7307.TF$M, file="../jun17.R.GSE7307.TF.M.tab", sep="\t")


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

# TODO: GENERATE BOOL DIST MATRIX
# load pre-computed boolean distance matrices
GSE31448.BOOL.D.row <- as.dist(as.matrix(read.table("~/brca/jun17.R.GSE31448.TF.BOOL.tab.booldist.tab", header=T, row.names=1, check.names=F)))
stopifnot(all(rownames(GSE31448.BOOL.D.row)==colnames(GSE31448.BOOL.D.row)))
stopifnot(all(rownames(GSE31448.BOOL.D.row)==rownames(GSE31448.TF$BOOL)))
GSE7307.BOOL.D.row <- as.dist(as.matrix(read.table("~/brca/jun17.R.GSE7307.TF.BOOL.tab.booldist.tab", header=T, row.names=1, check.names=F)))
stopifnot(all(rownames(GSE7307.BOOL.D.row)==colnames(GSE7307.BOOL.D.row)))
stopifnot(all(rownames(GSE7307.BOOL.D.row)==rownames(GSE7307.TF$BOOL)))

pdf("~/GSE7307.gsplom.tf.jun17.pdf", width=100, height=100)
R.GSE7307.TF <- splom(GSE7307.TF$BOOL, GSE7307.TF$DCOR, useRaster=T, sym=T, high.sat=T, MIN=0.2, MAX=1, row.cls.dist=GSE7307.BOOL.D.row)
dev.off()

pdf("~/GSE31448.gsplom.tf.jun17.pdf", width=100, height=100)
R.GSE31448.TF <- splom(GSE31448.TF$BOOL, GSE31448.TF$DCOR, useRaster=T, sym=T, high.sat=T, MIN=0.2, MAX=1, row.cls.dist=GSE31448.BOOL.D.row)
dev.off()

pdf("~/GSE7307.gsplom.tf.jun17.dendro.pdf", width=180, height=12)
plot(R.GSE7307.TF$Rhclust, main="GSE7307 TF+BRCA1/2")
dev.off()

pdf("~/GSE31448.gsplom.tf.jun17.dendro.pdf", width=180, height=12)
plot(R.GSE31448.TF$Rhclust, main="GSE31448 TF+BRCA1/2")
dev.off()

save(R.GSE7307.TF, R.GSE31448.TF, file="../jun17.BRCA.TF.gsploms.RData")

pdf("~/GSE31448.TF.summary.plots.jun17.pdf")
H <- summary.plots(GSE31448.TF$BOOL, GSE31448.TF$DCOR, sym=T)
dev.off()

pdf("~/GSE7307.TF.summary.plots.jun17.pdf")
H <- summary.plots(GSE7307.TF$BOOL, GSE7307.TF$DCOR, sym=T)
dev.off()

#
qq <- match(rownames(GSE31448.BOOL),rownames(GSE31448.DCOR))
pdf("~/GSE31448.summary.plots.jun17.pdf")
H <- summary.plots(GSE31448.BOOL[!is.na(qq),!is.na(qq)], GSE31448.DCOR[na.omit(qq),na.omit(qq)], sym=T)
dev.off()

qq <- match(rownames(GSE7307.BOOL),rownames(GSE7307.DCOR))
pdf("~/GSE7307.summary.plots.jun17.pdf")
H <- summary.plots(GSE7307.BOOL[!is.na(qq),!is.na(qq)], GSE7307.DCOR[na.omit(qq),na.omit(qq)], sym=T)
dev.off()
