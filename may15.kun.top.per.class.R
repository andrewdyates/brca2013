library(Biobase)
load("../STEPS.E7307.may14.R")
load("../BRCA1.dcor.RData")
load("../BRCA2.dcor.RData")
load("../BRCA.CLS.RData")
load("../E7307.may14.genelvl.exprfilt.RData")

qcm.gene.list <- as.character(read.table("../qcm_gene_list.txt")[,1])
# get adj matrix of PINA
# get transcription factor list
 
# read PINA list
# for BRCA1, BRCA2
# get list of all classes, ranked by decreasing dcor
# flag in qcm list, in transcription factor list, in PINA list
