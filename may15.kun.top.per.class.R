library(Biobase)
load("../STEPS.E7307.may14.R")
load("../BRCA1.dcor.RData")
load("../BRCA2.dcor.RData")
load("../BRCA.CLS.RData")
load("../E7307.may14.genelvl.exprfilt.RData")

load("~/tftargets/tf.adj.RData")
load("~/pina/pina.adj.RData")
source("~/pymod/dependency_glyph_splom/lib.R")
qcm.gene.list <- as.character(read.table("../qcm_gene_list.txt")[,1])

# read PINA list
# for BRCA1, BRCA2
# get list of all classes, ranked by decreasing dcor
# flag in qcm list, in transcription factor list, in PINA list
# PINA.ADJ, qcm.gene.list, TF.ADJ

all(rownames(BRCA1.CLS) == rownames(BRCA2.CLS)) # TRUE
all(rownames(BRCA1.DCOR) == rownames(BRCA2.DCOR)) # TRUE
all(rownames(BRCA1.CLS) == rownames(BRCA2.DCOR)) # TRUE
any(colnames(TF.ADJ) == "BRCA1") #FALSE (not tf)
any(colnames(TF.ADJ) == "BRCA2") #FALSE (not tf)

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
  R$ppi <- colnames(PINA.ADJ)[which(PINA.ADJ[rownames(PINA.ADJ) == name,]==1)]
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
  pina.g.i <- which(rownames(PINA.ADJ) == g)
  tf.g.i <- which(rownames(TF.ADJ) == g)
  dlist.sort <- sort(dlist, decreasing=T)
  # other gene is transcription factor
  is.trans <- function(s) s %in% colnames(TF.ADJ)
  # other gene is transcription factor, and it regulates this gene
  is.trans.reg <- function(s) s %in% colnames(TF.ADJ) && TF.ADJ[tf.g.i,which(s==colnames(TF.ADJ))]
  # other gene and this gene have ppi interaction
  is.ppi <- function(s) s %in% colnames(PINA.ADJ) && PINA.ADJ[pina.g.i,which(s==colnames(PINA.ADJ))]
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

BRCA1.R <- get.cls.lists(BRCA1.DCOR, BRCA1.CLS, "BRCA1")
BRCA2.R <- get.cls.lists(BRCA2.DCOR, BRCA2.CLS, "BRCA2")
save(BRCA1.R,BRCA2.R, file="../BRCA.R.may17.RData")

save.results(BRCA1.R, "BRCA1")
save.results(BRCA2.R, "BRCA2")
