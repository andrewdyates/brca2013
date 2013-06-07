#genes.e7307 <- readLines("../BRCA.network.gene.list.txt")
#genes.selected <- readLines("../selected_symbols_w_data.txt")
load("/home/ayates/brca/GSE31448/GSE31448.SCAN.pkl.DCOR.values.tab.RData")
DCOR.31448 <- M
load("/home/ayates/brca/GSE31448/GSE31448.may29.BOOL.RData")
CLS.31448 <- CLS
load("~/brca/GSE7307/jetset/E7307.may14.genelvl.exprfilt.nice.pkl.DCOR.tab.RData")
DCOR.7307 <- M
load("~/brca/GSE7307/jetset/GSE7307.may14.BOOL.RData")
CLS.7307 <- CLS
load("~/code/transcription_factors/all_tf_with_targs_adjm.tab.RData")
# load PPI
load("~/code/pina_parser/pina_compiled_adjm.tab.RData")
# load qcm-list
qcm.gene.list <- readLines("~/code/brca_qcm_list/qcm_gene_list_clean.txt")

source("/home/ayates/code/dependency_glyph_splom/lib.R")

pdf("all.genes.GSPLOM.summary.GSE31448.pdf")
summary.plots(CLS.31448, DCOR.31448, sym=T)
dev.off()

# check names
all(rownames(CLS.31448) == colnames(CLS.31448))
all(rownames(DCOR.31448) == colnames(DCOR.31448)) #F
# remove "X" in front of column names
colnames(DCOR.31448) <- sub("X","",colnames(DCOR.31448))
all(rownames(DCOR.31448) == colnames(DCOR.31448)) #T
all(colnames(DCOR.31448) == colnames(CLS.31448))
all(rownames(DCOR.31448) == rownames(CLS.31448))

all(rownames(CLS.31448) == rownames(DCOR.31448))
all(colnames(CLS.31448) == colnames(DCOR.31448))
all(rownames(CLS.7307) == colnames(CLS.7307))
all(rownames(CLS.7307) == rownames(DCOR.7307))
all(colnames(CLS.7307) == colnames(DCOR.7307)) # F
# somewhere, R converted "." to "-" in column names for DCOR. Fix this.
colnames(DCOR.7307) <- colnames(CLS.7307)

# select BRCA1, BRCA2, plot summaries
brca1.i <- which(rownames(CLS.31448)=='672') #14449
brca2.i <- which(rownames(CLS.31448)=='675') #14475

# fix object type
DCOR.31448 <- as.matrix(DCOR.31448)
DCOR.7307 <- as.matrix(DCOR.7307)

pdf("BRCA1.GSPLOM.summary.GSE31448.pdf")
summary.plots.vector(CLS.31448[brca1.i,], DCOR.31448[brca1.i,])
dev.off()
pdf("BRCA2.GSPLOM.summary.GSE31448.pdf")
summary.plots.vector(CLS.31448[brca2.i,], DCOR.31448[brca2.i,])
dev.off()

# convert entrezg to gene symbol in 31448
library(hgu133plus2hsentrezg.db)
IDS <- paste(rownames(CLS.31448),"_at", sep="")
SYMS <- mget(IDS, hgu133plus2hsentrezgSYMBOL, ifnotfound=NA)
rownames(CLS.31448) <- as.character(SYMS)
colnames(CLS.31448) <- as.character(SYMS)
rownames(DCOR.31448) <- as.character(SYMS)
colnames(DCOR.31448) <- as.character(SYMS)



get.cls.lists <- function(DCOR, CLS, name) {
  R <- list()
  Z <- split(DCOR, CLS)
  R$table <- lapply(Z, function(z) summarize.list(z,name))
                                        # TF Targets (to name)
  R$tft <- colnames(TF.ADJ)[which(TF.ADJ[rownames(TF.ADJ) == name,]==1)]
  R$tft.dcor <- DCOR[which(names(DCOR) %in% R$tft)]
  R$tft.cls <- CLS[which(names(CLS) %in% R$tft)]
  R$tft.miss <- R$tft[which(!R$tft %in% names(DCOR))]
                                        # All TF
  R$tf <- colnames(TF.ADJ)
  R$tf.dcor <- DCOR[which(names(DCOR) %in% R$tf)]
  R$tf.cls <- CLS[which(names(CLS) %in% R$tf)]
  R$tf.miss <- R$tf[which(!R$tf %in% names(DCOR))]
                                        # PINA PPI
  R$ppi <- colnames(PPI.ADJ)[which(PPI.ADJ[rownames(PPI.ADJ) == name,]==1)]
  R$ppi.dcor <- DCOR[which(names(DCOR) %in% R$ppi)]
  R$ppi.cls <- CLS[which(names(CLS) %in% R$ppi)]
  R$ppi.miss <- R$ppi[which(!R$ppi %in% names(DCOR))]
                                        # QCM cluster lists
  R$qcm <- qcm.gene.list
  R$qcm.dcor <- DCOR[names(DCOR) %in% qcm.gene.list]
  R$qcm.cls <- CLS[names(CLS) %in% qcm.gene.list]
  R$qcm.miss <- R$qcm[!R$qcm %in% names(DCOR)]
  R
}

summarize.list <- function(dlist, g="BRCA1") {
  # we already know that BRCA1 and BRCA2 are in the rows of both PINA and TF ADJ matrices
  R <- list()
  pina.g.i <- which(rownames(PPI.ADJ) == g)
  tf.g.i <- which(rownames(TF.ADJ) == g)
  dlist.sort <- sort(dlist, decreasing=T)
  # other gene is transcription factor
  is.trans <- function(s) s %in% colnames(TF.ADJ)
  # other gene is transcription factor, and it regulates this gene
  is.trans.reg <- function(s) s %in% colnames(TF.ADJ) && TF.ADJ[tf.g.i,which(s==colnames(TF.ADJ))]
  # other gene and this gene have ppi interaction
  is.ppi <- function(s) s %in% colnames(PPI.ADJ) && PPI.ADJ[pina.g.i,which(s==colnames(PPI.ADJ))]
  # is in QCM list
  is.qcm <- function(s) s %in% qcm.gene.list

  trans.v <- sapply(names(dlist.sort), is.trans)
  trans.targ.v <- sapply(names(dlist.sort), is.trans.reg)
  ppi.v <- sapply(names(dlist.sort), is.ppi)
  qcm.v <- sapply(names(dlist.sort), is.qcm)
  data.frame(dlist.sort, trans.v, trans.targ.v, ppi.v, qcm.v)
}

save.results <- function(R, name) {
  for (cls in names(R$table)) {
    fname <- paste0(name,".",cls,".tab")
    write.table(R$table[[cls]], file=fname, sep="\t")
  }
  R
  plot.summary(R$tf.dcor, R$tf.cls, R$tf.miss, paste0(name,".tf.summ.pdf"))
  plot.summary(R$tft.dcor, R$tft.cls, R$tft.miss, paste0(name,".tft.summ.pdf"))
  plot.summary(R$ppi.dcor, R$ppi.cls, R$ppi.miss, paste0(name,".ppi.summ.pdf"))
  plot.summary(R$qcm.dcor, R$qcm.cls, R$qcm.miss, paste0(name,".qcm.summ.pdf"))
}

plot.summary <- function(DCOR, CLS, MISS, fname) {
  DCOR <- c(DCOR, rep(0,length(MISS)))
  CLS <- c(CLS, rep(0,length(MISS)))
  pdf(fname)
  summary.plots.vector(CLS, DCOR)
  dev.off()
}

BRCA1.31448.R <- get.cls.lists(DCOR.31448[brca1.i,], CLS.31448[brca1.i,], "BRCA1")
BRCA2.31448.R <- get.cls.lists(DCOR.31448[brca2.i,], CLS.31448[brca2.i,], "BRCA2")

#CLS.7307
brca1.i2 <- which(rownames(DCOR.7307)=="BRCA1")
brca2.i2 <- which(rownames(DCOR.7307)=="BRCA2")
BRCA1.7307.R <- get.cls.lists(DCOR.7307[brca1.i2,], CLS.7307[brca1.i2,], "BRCA1")
BRCA2.7307.R <- get.cls.lists(DCOR.7307[brca2.i2,], CLS.7307[brca2.i2,], "BRCA2")

save(BRCA1.31448.R, BRCA2.31448.R, BRCA1.7307.R, BRCA2.7307.R, file="../may31.7307.31448.brca.cls.lists.RData")


save.results(BRCA1.31448.R, "BRCA1.31448")
save.results(BRCA2.31448.R, "BRCA2.31448")
save.results(BRCA1.7307.R, "BRCA1.7307")
save.results(BRCA2.7307.R, "BRCA2.7307")
