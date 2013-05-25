load("~/brca/GSE7307.may14.BOOL.RData")
REMAP <- read.table("E7307.genesym.remap.txt", stringsAsFactors=F)[,2]
rownames(CLS) <- REMAP
colnames(CLS) <- REMAP
CLS.row.D <- gen.glyph.dist.m(CLS)
load("~/brca/GSE7307.may14.BOOL.row.D.RData")
