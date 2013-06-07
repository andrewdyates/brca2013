BRCA1.31448.1 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA1.31448.1.tab")
BRCA1.31448.2 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA1.31448.2.tab")
BRCA1.31448.3 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA1.31448.3.tab")
BRCA1.31448.4 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA1.31448.4.tab")
BRCA2.31448.1 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA2.31448.1.tab")
BRCA2.31448.2 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA2.31448.2.tab")
BRCA2.31448.3 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA2.31448.3.tab")
BRCA2.31448.4 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA2.31448.4.tab")
BRCA1.7307.1 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA1.7307.1.tab")
BRCA1.7307.2 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA1.7307.2.tab")
BRCA1.7307.3 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA1.7307.3.tab")
BRCA1.7307.4 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA1.7307.4.tab")
BRCA2.7307.1 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA2.7307.1.tab")
BRCA2.7307.2 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA2.7307.2.tab")
BRCA2.7307.3 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA2.7307.3.tab")
BRCA2.7307.4 <- read.table("/Users/z/Dropbox/biostat/brca/may31.brca.7307.31448.plots/BRCA2.7307.4.tab")

# BRCA1
brca1.1 <- intersect(rownames(BRCA1.31448.1), rownames(BRCA1.7307.1))
brca1.1
brca1.1.tab <- BRCA1.31448.1[rownames(BRCA1.31448.1) %in% brca1.1,]
BRCA1.7307.1[rownames(BRCA1.7307.1) %in% brca1.1,]
write.table(brca1.1.tab, file="../may31.mix+bc.brca1.1.tab", sep="\t", col.names=NA)

brca1.2 <- intersect(rownames(BRCA1.31448.2), rownames(BRCA1.7307.2)) #none
brca1.2

brca1.3 <- intersect(rownames(BRCA1.31448.3), rownames(BRCA1.7307.3))
brca1.3
brca1.3.tab <- BRCA1.31448.3[rownames(BRCA1.31448.3) %in% brca1.3,]
BRCA1.7307.3[rownames(BRCA1.7307.3) %in% brca1.3,]
write.table(brca1.3.tab, file="../may31.mix+bc.brca1.3.tab", sep="\t", col.names=NA)

#brca1.4 <- intersect(rownames(BRCA1.31448.4), rownames(BRCA1.7307.4))
#brca1.4.tab <- BRCA1.7307.4[rownames(BRCA1.7307.4) %in% brca1.4,]

# BRCA2
brca2.1 <- intersect(rownames(BRCA2.31448.1), rownames(BRCA2.7307.1))
brca2.1
brca2.1.tab <- BRCA2.31448.1[rownames(BRCA2.31448.1) %in% brca2.1,]
BRCA2.7307.1[rownames(BRCA2.7307.1) %in% brca2.1,]
write.table(brca2.1.tab, file="../may31.mix+bc.brca2.1.tab", sep="\t", col.names=NA)

brca2.2 <- intersect(rownames(BRCA2.31448.2), rownames(BRCA2.7307.2)) #none
brca2.2

brca2.3 <- intersect(rownames(BRCA2.31448.3), rownames(BRCA2.7307.3))
length(brca2.3) # 74
brca2.3.tab <- BRCA2.31448.3[rownames(BRCA2.31448.3) %in% brca2.3,]
BRCA2.7307.3[rownames(BRCA2.7307.3) %in% brca2.3,]
write.table(brca2.3.tab, file="../may31.mix+bc.brca2.3.tab", sep="\t", col.names=NA)

brca12.1 <- intersect(brca1.1, brca2.1)
brca12.1.tab <- BRCA1.31448.1[rownames(BRCA1.31448.1) %in% brca12.1,]
write.table(brca2.1.tab, file="../may31.mix+bc.brca12.1.tab", sep="\t", col.names=NA)

brca12.3 <- intersect(brca1.3, brca2.3)
brca12.3.tab <- BRCA1.31448.3[rownames(BRCA1.31448.3) %in% brca12.3,]
write.table(brca12.3.tab, file="../may31.mix+bc.brca12.3.tab", sep="\t", col.names=NA)

