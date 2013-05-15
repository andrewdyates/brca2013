library(Biobase)
library(jetset)
library(org.Hs.eg.db)
load("../E7307.may14.RData")

## TEST jetset module
## ----------------------------------------
length(na.omit(unique(scores.hgu133plus2$EntrezID)))
jmap('hgu133plus2', symbol='BRCA2')
#      BRCA2 
#"214727_at" 
jscores('hgu133plus2', symbol='BRCA2')
##             nProbes EntrezID process specificity coverage    robust   overall
## 214727_at        11      675     296   0.9090909        1 0.6103365 0.5548513
## 208368_s_at      11      675     836   0.9090909        1 0.2479580 0.2254164
##             symbol
## 214727_at    BRCA2
## 208368_s_at  BRCA2
all.jscores <- jscores('hgu133plus2')
all.entrezID <- unique(na.omit(all.jscores$EntrezID))
all.entrezID.jmap <- jmap('hgu133plus2', eg=all.entrezID)
all.entrezID.genes <- mget(names(all.entrezID.jmap), org.Hs.egSYMBOL, ifnotfound=NA)

qq <- match(all.entrezID.jmap, rownames(exprs(E7307)))
E7307.genes <- E7307[qq,]
all(rownames(exprs(E7307.genes)) == all.entrezID.jmap) #TRUE
all(names(all.entrezID.jmap) == names(all.entrezID.genes)) #TRUE

featureData(E7307.genes)$jetset.entrezID <- names(all.entrezID.jmap)
featureData(E7307.genes)$jetset.genesym <- all.entrezID.genes
sum(is.na(featureData(E7307.genes)$jetset.genesym))
[1] 49

save(E7307.genes, file="../E7307.may14.genelvl.RData")
