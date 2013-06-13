load("../jun12.BRCA.dep.matrices.TF.RData")
source("~/code/dependency_glyph_splom/lib.R")
# write matrices to file
write.table(GSE31448.TF$BOOL, file="../jun12.R.GSE31448.TF.BOOL.tab", sep="\t", quote=F)
write.table(GSE31448.TF$WEAK, file="../jun12.R.GSE31448.TF.WEAK.tab", sep="\t", quote=F)
write.table(GSE31448.TF$DCOR, file="../jun12.R.GSE31448.TF.DCOR.tab", sep="\t", quote=F)

write.table(GSE7307.TF$BOOL, file="../jun12.R.GSE7307.TF.BOOL.tab", sep="\t", quote=F)
write.table(GSE7307.TF$WEAK, file="../jun12.R.GSE7307.TF.WEAK.tab", sep="\t", quote=F)
write.table(GSE7307.TF$DCOR, file="../jun12.R.GSE7307.TF.DCOR.tab", sep="\t", quote=F)

# ==================================================
# Transcription Factor Network
# ==================================================
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

# load pre-computed boolean distance matrices
GSE31448.BOOL.D.row <- as.dist(as.matrix(read.table("~/brca/jun12.R.GSE31448.TF.BOOL.tab.booldist.tab", header=T, row.names=1, check.names=F)))
stopifnot(all(rownames(GSE31448.BOOL.D.row)==colnames(GSE31448.BOOL.D.row)))
stopifnot(all(rownames(GSE31448.BOOL.D.row)==rownames(GSE31448.TF$BOOL)))
GSE7307.BOOL.D.row <- as.dist(as.matrix(read.table("~/brca/jun12.R.GSE7307.TF.BOOL.tab.booldist.tab", header=T, row.names=1, check.names=F)))
stopifnot(all(rownames(GSE7307.BOOL.D.row)==colnames(GSE7307.BOOL.D.row)))
stopifnot(all(rownames(GSE7307.BOOL.D.row)==rownames(GSE7307.TF$BOOL)))



pdf("~/GSE7307.gsplom.tf.pdf", width=100, height=100)
R.GSE7307.TF <- splom(GSE7307.TF$BOOL, GSE7307.TF$DCOR, useRaster=T, sym=T, high.sat=T, MIN=0.2, MAX=1, row.cls.dist=GSE7307.BOOL.D.row)
dev.off()

pdf("~/GSE31448.gsplom.tf.pdf", width=100, height=100)
R.GSE31448.TF <- splom(GSE31448.TF$BOOL, GSE31448.TF$DCOR, useRaster=T, sym=T, high.sat=T, MIN=0.2, MAX=1, row.cls.dist=GSE31448.BOOL.D.row)
dev.off()

pdf("~/GSE31448.gsplom.tf.dcorweight0.1.pdf", width=100, height=100)
R.GSE31448.TF.dcor0.1 <- splom(GSE31448.TF$BOOL, GSE31448.TF$DCOR, useRaster=T, sym=T, high.sat=T, MIN=0.2, MAX=1, row.cls.dist=GSE31448.BOOL.D.row,DCOR.weight=0.1)
dev.off()

pdf("~/GSE31448.gsplom.tf.dcorweight40.pdf", width=100, height=100)
R.GSE31448.TF.dcor40 <- splom(GSE31448.TF$BOOL, GSE31448.TF$DCOR, useRaster=T, sym=T, high.sat=T, MIN=0.2, MAX=1, row.cls.dist=GSE31448.BOOL.D.row,DCOR.weight=40)
dev.off()

pdf("~/GSE7307.gsplom.tf.dendro.pdf", width=180, height=12)
plot(R.GSE7307.TF$Rhclust, main="GSE7307 TF+BRCA1/2")
dev.off()
pdf("~/GSE31448.gsplom.tf.dendro.pdf", width=180, height=12)
plot(R.GSE31448.TF$Rhclust, main="GSE31448 TF+BRCA1/2")
dev.off()


save(R.GSE7307.TF, R.GSE31448.TF, file="../jun12.BRCA.TF.gsploms.RData")
