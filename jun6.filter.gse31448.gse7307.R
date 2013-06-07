#bucki paths
library(Biobase)
library(hgu133plus2hsentrezg.db)

load("~/brca/GSE31448/normed/GSE31448.SCANUPC.RData")
load("~/brca/GSE7307/normed/GSE7307.SCANUPC.RData")
AttrT.31448 <- as.data.frame(t(read.table("~/brca/GSE31448/geo_download/GSE31448_GPL570.samples.tab", sep="\t", row.names=1, header=T, stringsAsFactors=F)))
AttrT.7307 <- as.data.frame(t(read.table("~/brca/GSE7307/raw/GSE7307_GPL570.samples.tab", sep="\t", row.names=1, header=T, stringsAsFactors=F)))

names(AttrT.31448)
##  [1] "erbb2"              "label_protocol"     "pr.ihc"            
##  [4] "mfs"                "type"               "relation"          
##  [7] "pr ihc status"      "source_name"        "title"             
## [10] "scan_protocol"      "data_processing"    "supplementary_file"
## [13] "prihc"              ".tab ihc status"    "description"       
## [16] "pt"                 "sbr grade"          "histology"         
## [19] " ihc status"        "egfr ihc status"    "age.at.diagnosis"  
## [22] "mfsdel (month)"     "geo_accession"      "pn"                
## [25] "er ihc status"      "p53 ihc"            "er.ihc"            
## [28] "dfs evt"            "top2a ihc status"   "dfs time (months)" 
## [31] "ki67 ihc"           "foxa1 ihc status"   ".ht ihc status"    
## [34] "grade sbr"          "igf1r ihc status"   "p ihc status"      
## [37] "hyb_protocol"       "molecular subtype" 
names(AttrT.7307)
 ## [1] "relation"                        "title"                          
 ## [3] "Disease/Normal or Treatment [C]" "description"                    
 ## [5] "supplementary_file"              "Disease type"                   
 ## [7] "last_update_date"                "Gender"                         
 ## [9] "geo_accession"                   "Tissue/Cell Line [C]"           

dim(S.GSE31448)
#[1] 18896   357
dim(S.GSE7307)
#[1] 18896   677

# get b
stds.31448 <- apply(S.GSE31448,1,sd)
stds.7307 <- apply(S.GSE7307,1,sd)
b.31448 <- as.numeric(quantile(stds.31448, 0.03)*2)
#[1] 0.09676083
b.7307 <- as.numeric(quantile(stds.7307, 0.03)*2)
#[1] 0.1541373
sum(stds.31448 > b.31448)
#[1] 14794
sum(stds.7307 > b.7307)
#[1] 13766

# filter by variance
S.GSE31448 <- S.GSE31448[stds.31448 > b.31448,]
S.GSE7307 <- S.GSE7307[stds.7307 > b.7307,]

# save results
save(S.GSE31448, S.GSE7307, b.31448, b.7307, file="../jun6.GSE7307.GSE31448.select.RData")
write.table(S.GSE31448, file="../jun6.GSE31448.select.tab", sep="\t")
write.table(S.GSE7307, file="../jun6.GSE7307.select.tab", sep="\t")
