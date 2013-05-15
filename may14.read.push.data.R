library("Biobase")
library(jetset)
gse7307.norm <- read.table(file="../7307_scanupc_sharpnack.txt", sep="\t", header=T, row.names=1, comment="", quote="")

#E7307 <- ExpressionSet(as.matrix(gse7307.norm))
#save(E7307, file="../gse7307.SCAN.may14.RData")
load(file="../gse7307.SCAN.may14.RData")
GPL570 <- read.table(file="../GPL570-13270.txt", sep="\t", quote="", header=T, row.names=1, as.is=T)
qq <- match(rownames(exprs(E7307)), rownames(GPL570))
all(qq == 1:length(qq)) # FALSE!
pData(featureData(E7307)) <- GPL570[qq,]
#all(match(rownames(exprs(E7307)), rownames(GPL570)), 1:dim(GPL570)[1]) # TRUE.
varLabels(featureData(E7307))
##  [1] "GB_ACC"                           "SPOT_ID"                         
##  [3] "Species.Scientific.Name"          "Annotation.Date"                 
##  [5] "Sequence.Type"                    "Sequence.Source"                 
##  [7] "Target.Description"               "Representative.Public.ID"        
##  [9] "Gene.Title"                       "Gene.Symbol"                     
## [11] "ENTREZ_GENE_ID"                   "RefSeq.Transcript.ID"            
## [13] "Gene.Ontology.Biological.Process" "Gene.Ontology.Cellular.Component"
## [15] "Gene.Ontology.Molecular.Function"

save(E7307, file="../E7307.may14.RData")
