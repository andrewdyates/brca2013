library(Biobase)
load("../E7307.may14.genelvl.exprfilt.RData")
writeLines(rownames(M), "E7307.original.row.names.txt")
