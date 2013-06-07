#install.packages("hgu133plus2hsentrezgprobe_17.1.0.tar.gz", repos=NULL, type="source")
#install.packages("hgu133plus2hsentrezg.db_17.1.0.tar.gz", repos=NULL, type="source")
library(hgu133plus2hsentrezgprobe)
library(hgu133plus2hsentrezg.db)
library(SCAN.UPC)
celFilePath = "GSM887367.CEL.gz"
norm = SCAN(celFilePath)
norm.cust = SCAN(celFilePath, probeSummaryPackage=hgu133plus2hsentrezgprobe)

IDS <- rownames(exprs(norm.cust)) # note: IDS are entrez_ids appended by _at
EIDS <- mget(IDS, hgu133plus2hsentrezgENTREZID, ifnotfound=NA)
rm.noann <- is.na(EIDS)
norm.filt <- norm.cust[!rm.noann,]
rownames(norm.filt) <- EIDS[!rm.noann]
# every row in EID is now uniquely identified by an entrez ID
