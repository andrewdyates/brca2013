library(Biobase)
load("~/brca/E7307.may14.genelvl.exprfilt.nice.pkl.DCOR.tab.RData")
DCOR <- M
load("~/brca/GSE7307.may14.BOOL.RData")
load("~/tftargets/tf.adj.may22.2013.RData")
load("~/pina/pina.may22.2013.adj.RData")
load("../E7307.may14.genelvl.exprfilt.RData")
source("~/pymod/dependency_glyph_splom/lib.R")
qcm <- readLines("~/pymod/brca_qcm_list/qcm_gene_list_clean.txt")
select.syms <- readLines("BRCA.network.gene.list.txt")
sym.remap <- read.table("E7307.genesym.remap.txt", stringsAsFactors=F)

# first pass summary plots
# NOTE: too many outliers drawn distorts per-class boxplots
pdf("gsplom.summary.all.pdf")
summary.plots(CLS, DCOR, sym=T)
dev.off()


# TODO
# Load weaks
load("../E7307.may22.genelvl.exprfilt.WEAK.RData")
# Load PCC
#load(?)


# Verify ID mappings

# summary plots for all pairs DCOR, CLS

# this fails for a large number of outliers?

all(rownames(CLS)==rownames(DCOR))


