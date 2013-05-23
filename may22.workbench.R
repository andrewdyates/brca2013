library(Biobase)
load("../E7307.may14.genelvl.exprfilt.RData")
source("~/pymod/dependency_glyph_splom/lib.R")

REMAP <- read.table("E7307.genesym.remap.txt", stringsAsFactors=F)
featureData(E)$bestsym <- REMAP[,2]
rownames(M) <- REMAP[,2]

save(M,E,b,file="../E7307.may22.genelvl.exprfilt.RData")
