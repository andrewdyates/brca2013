source("~/pymod/dependency_glyph_splom/lib.R")
library(Biobase)
load("../GSE7307.select.may25.RData")
load("~/tftargets/tf.adj.may22.2013.RData")
load("~/pina/pina.may22.2013.adj.RData")

pdf("select.summary.plots.pdf")
summary.plots(CLS.s, DCOR.s, sym=T)
dev.off()
summary(as.factor(CLS.s))
#      0       1       2       3       4       5       6       7 
#     16  123409   16291  123409 1281772  269422       6      84

Z7 <- which(CLS.s==7, arr.ind=T)
Z6 <- which(CLS.s==6, arr.ind=T)

rownames(CLS.s)[Z6[,1]]
rownames(CLS.s)[Z6[,2]]

pdf("antigsplom.pdf")
splom(CLS.s[Z6[,1],Z6[,2]], DCOR.s[Z6[,1],Z6[,2]], reorder=F, MAX=1, MIN=0.15, asGlyphs=T, lwd=3)
splom(CLS.s[Z6[,1],Z6[,2]], DCOR.s[Z6[,1],Z6[,2]], reorder=F, MAX=1, MIN=0.15, asGlyphs=F, grid=F)
splom(CLS.s[Z6[,1],Z6[,2]], DCOR.s[Z6[,1],Z6[,2]], reorder=T, MAX=1, MIN=0.15, asGlyphs=F, grid=F)
splom(CLS.s[Z6[,1],Z6[,2]], DCOR.s[Z6[,1],Z6[,2]], reorder=T, MAX=1, MIN=0.15, asGlyphs=T, grid=F, high.sat=F)
splom(CLS.s[Z6[,1],Z6[,2]], DCOR.s[Z6[,1],Z6[,2]], reorder=T, MAX=1, MIN=0.15, asGlyphs=T, lwd=3, high.sat=F)
dev.off()

CLS.s[Z6[,1],Z6[,2]]
#[1] "SOX2"  "REST"  "SIRT3" "RREB1" "SOX8"  "YBX1" 
colnames(CLS.s)[Z6[,2]]
#[1] "REST"  "SOX2"  "YBX1"  "SOX8"  "RREB1" "SIRT3"

pdf("Selected.GSPLOM.pdf", width=100, height=100)
R <- splom(CLS.s, DCOR.s, grid=F, sym=T, asGlyphs=F, useRaster=T)
dev.off()
save(R, file="../Selected.GSPLOM.Rdata")
