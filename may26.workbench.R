source("~/pymod/dependency_glyph_splom/lib.R")
library(Biobase)
load("../Selected.GSPLOM.Rdata")
load("../GSE7307.select.may25.RData")

pdf("Selected.GSPLOM.dendro.pdf", width=200, height=10)
plot(R$Rhclust)
dev.off()

Rhclust <- as.dendrogram(hclust(R$D.row, method="average"))
Rowv <- rowMeans(DCOR.s, na.rm = TRUE)
Rhclust <- reorder(Rhclust, Rowv)
qq <- order.dendrogram(Rhclust)


pdf("Selected.GSPLOM.avg.pdf", width=100, height=100)
R.avg <- splom(CLS.s[qq,qq], DCOR.s[qq,qq], grid=F, sym=T, asGlyphs=F, useRaster=T, reorder=F)
dev.off()

pdf("Selected.GSPLOM.dendro.avg.pdf", width=200, height=10)
plot(Rhclust)
dev.off()

# compare coherence curves for avg, complete linkage
H.complete <- as.hclust(R$Rhclust)
H.avg <- as.hclust(Rhclust)
# This takes forever due to naive implementation
#Z.complete <- get.compression(CLS.s,H.complete,as.matrix(DCOR.s),min.dcor=0.2)

## C.complete.100 <- collapse.cls(CLS.s,cutree(H.complete,100),as.matrix(DCOR.s))
## S.complete.100 <- get.coh.M.score(C.complete.100, 0.2)
## C.avg.100 <- collapse.cls(CLS.s,cutree(H.avg,100),as.matrix(DCOR.s))
## S.avg.100 <- get.coh.M.score(C.avg.100, 0.2)
## S.complete.100$edge.flaws / S.complete.100$edge.all.n #[1] 0.1525253
## S.avg.100$edge.flaws / S.avg.100$edge.all.n           #[1] 0.0979798

## C.complete.150 <- collapse.cls(CLS.s,cutree(H.complete,150),as.matrix(DCOR.s))
## S.complete.150 <- get.coh.M.score(C.complete.150, 0.2)
## C.avg.150 <- collapse.cls(CLS.s,cutree(H.avg,150),as.matrix(DCOR.s))
## S.avg.150 <- get.coh.M.score(C.avg.150, 0.2)
## S.complete.150$edge.flaws / S.complete.150$edge.all.n #[1] 0.08170022
## S.avg.150$edge.flaws / S.avg.150$edge.all.n           #[1] 0.05852349

## C.complete.200 <- collapse.cls(CLS.s,cutree(H.complete,200),as.matrix(DCOR.s))
## S.complete.200 <- get.coh.M.score(C.complete.200, 0.2)
## C.avg.200 <- collapse.cls(CLS.s,cutree(H.avg,200),as.matrix(DCOR.s))
## S.avg.200 <- get.coh.M.score(C.avg.200, 0.2)
## S.complete.200$edge.flaws / S.complete.200$edge.all.n #[1] 0.04839196
## S.avg.200$edge.flaws / S.avg.200$edge.all.n           #[1] 0.03562814
## S.avg.200$edge.flaws / S.avg.200$edge.n               #[1] 0.05981608

## C.complete.250 <- collapse.cls(CLS.s,cutree(H.complete,250),as.matrix(DCOR.s))
## S.complete.250 <- get.coh.M.score(C.complete.250, 0.2)
## C.avg.250 <- collapse.cls(CLS.s,cutree(H.avg,250),as.matrix(DCOR.s))
## S.avg.250 <- get.coh.M.score(C.avg.250, 0.2)
## S.complete.250$edge.flaws / S.complete.250$edge.all.n #[1] 0.03261044
## S.avg.250$edge.flaws / S.avg.250$edge.all.n           #[1] 0.02682731
## S.complete.250$edge.flaws / S.complete.250$edge.n     #[1] 0.05056796
## S.avg.250$edge.flaws / S.avg.250$edge.n               #[1] 0.04427594
## S.avg.250$edge.n / S.avg.250$edge.all.n               #[1] 0.6059116

## C.avg.225 <- collapse.cls(CLS.s,cutree(H.avg,225),as.matrix(DCOR.s))
## S.avg.225 <- get.coh.M.score(C.avg.225, 0.2)
## S.avg.225$edge.flaws / S.avg.225$edge.n
## #[1] 0.05191221

## C.avg.230 <- collapse.cls(CLS.s,cutree(H.avg,230),as.matrix(DCOR.s))
## S.avg.230 <- get.coh.M.score(C.avg.230, 0.2)
## S.avg.230$edge.flaws / S.avg.230$edge.n
## #[1] 0.05106437

# AVERAGE LINKAGE, 5% extant, K=235
C.avg.235 <- collapse.cls(CLS.s,cutree(H.avg,235),as.matrix(DCOR.s))
S.avg.235 <- get.coh.M.score(C.avg.235, 0.2)
S.avg.235$edge.flaws / S.avg.235$edge.n
#[1] 0.04943588
S.avg.235$edge.flaws / S.avg.235$edge.all.n
#[1] 0.02964175
S.avg.235$edge.n / S.avg.235$edge.all.n
#[1] 0.5995999


## C.avg.233 <- collapse.cls(CLS.s,cutree(H.avg,233),as.matrix(DCOR.s))
## S.avg.233 <- get.coh.M.score(C.avg.233, 0.2)
## S.avg.233$edge.flaws / S.avg.233$edge.n
## #[1] 0.0504689
## C.avg.232 <- collapse.cls(CLS.s,cutree(H.avg,232),as.matrix(DCOR.s))
## S.avg.232 <- get.coh.M.score(C.avg.232, 0.2)
## S.avg.232$edge.flaws / S.avg.232$edge.n
## #[1] 0.05013029

## C.avg.300 <- collapse.cls(CLS.s,cutree(H.avg,300),as.matrix(DCOR.s))
## S.avg.300 <- get.coh.M.score(C.avg.300, 0.2)          # [1] 0.03070961
## S.avg.300$edge.flaws / S.avg.300$edge.all.n           # [1] 0.0187291
## S.avg.300$edge.n / S.avg.300$edge.all.n               # [1] 0.6098774
# ------------------------------

# HOW TO SAVE MEMBER LIST?
# TODO: all collapsed weaks
# 0: no class; 1: and; 2: rn4c (row necessary for col); 3: cn4r (col necessary for row); 4: xor; 5: mix, 6: no class
W.avg.235 <- collapse.weak(WEAK.s,cutree(H.avg,235))
summary(as.factor(WEAK.s))
#      1       2       3       4       5 
#  44873  188435  188435   31052 1361614 


save(C.avg.235, W.avg.235, file="../may26.k235.avglink.5pct.extant.collapsed.RData")
write.table(C.avg.235$CLS, file="../may26.k235.avglink.5pct.extant.collapsed.CLS.tab", sep="\t", col.names=NA)
write.table(C.avg.235$DCOR, file="../may26.k235.avglink.5pct.extant.collapsed.DCOR.tab", sep="\t", col.names=NA)
write.table(W.avg.235, file="../may26.k235.avglink.5pct.extant.collapsed.WEAK.tab", sep="\t", col.names=NA)
write.table(cutree(H.avg,235), file="../may26.k235.avglink.5ct.extant.collapsed.hclust.csv", sep=",", col.names=NA)

pdf("../may26.k235.avglink.5pct.extant.collapsed.GSPLOM.pdf", width=30, height=30)
R.k235 <- splom(C.avg.235$CLS, C.avg.235$DCOR, grid=F, sym=T)
dev.off()
pdf("../may26.k235.avglink.5pct.extant.collapsed.GSPLOM.dendro.pdf", width=30, height=30)
plot(R.k235$Rhclust)
dev.off()
