library(Biobase)
load("../E7307.may14.genelvl.exprfilt.RData")
source("~/pymod/dependency_glyph_splom/lib.R")

# remap gene symbols in expression matrix
REMAP <- read.table("E7307.genesym.remap.txt", stringsAsFactors=F)
featureData(E)$bestsym <- REMAP[,2]
rownames(M) <- REMAP[,2]
save(M,E,b,file="../E7307.may22.genelvl.exprfilt.RData")

load("~/tftargets/tf.adj.may22.2013.RData")
load("~/pina/pina.may22.2013.adj.RData")
# generate BRCA1+2 list for network
qcm <- as.character(read.table("~/pymod/brca_qcm_list/qcm_gene_list_clean.txt", stringsAsFactors=F)[,1]) #412
tfs <- colnames(TF.ADJ) # [1] 1375 (only "high quality", + lit)

brca1.i <- which(rownames(PINA.ADJ)=="BRCA1")
brca2.i <- which(rownames(PINA.ADJ)=="BRCA2")
ppi.1 <- names(which(PINA.ADJ[brca1.i,]==1)) # 132
ppi.2 <- names(which(PINA.ADJ[brca2.i,]==1)) # 56

net.list <- c(qcm, tfs, ppi.1, ppi.2, c("BRCA1","BRCA2"))
length(net.list) #1977
length(unique(net.list)) #1891
net.list <- sort(unique(net.list))
writeLines(net.list, "BRCA.network.gene.list.txt")
