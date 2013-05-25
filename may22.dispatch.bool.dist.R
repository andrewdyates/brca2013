load("~/brca/GSE7307.may14.BOOL.RData")
source("~/pymod/dependency_glyph_splom/lib.R")
REMAP <- read.table("~/brca/brca2013/E7307.genesym.remap.txt", stringsAsFactors=F)[,2]
rownames(CLS) <- REMAP
colnames(CLS) <- REMAP
CLS.row.D <- gen.glyph.dist.m(CLS)
save(CLS.row.D, CLS, file="~/brca/GSE7307.may14.BOOL.row.D.RData")
