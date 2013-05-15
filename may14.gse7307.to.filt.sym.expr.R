library(Biobase)
load("../E7307.may14.genelvl.RData")

# filter all probes without gene symbol
qq <- !is.na(featureData(E7307.genes)$jetset.genesym)
E <- E7307.genes[qq,]
# filter all probes without at least 4 probes expressed above 0.2
num.expr <- apply(exprs(E)>0.2, 1, sum)

qq <- which(featureData(E)$jetset.genesym %in% c("BRCA1", "BRCA2"))
num.expr[qq]
## 214727_at 204531_s_at @ 0.5
##         8          57
##  214727_at 204531_s_at @ 0.2
##         34         131 
std.expr <- apply(exprs(E), 1, sd)
std.expr.ranks <- rank(std.expr)
std.expr.ranks[qq] / length(std.expr)
## 214727_at 204531_s_at 
## 0.2844372   0.3997073
std.expr[qq]
#  214727_at 204531_s_at 
#  0.1657502   0.2148387 

b <- quantile(std.expr, 0.03)
#0.07683501 
sum(std.expr < b*2)
## [1] 4904
sum(std.expr < b*2)/length(std.expr)
#[1] 0.2563647 OK, use 25% variance cut
sum(num.expr == 0)
## [1] 939
sum(std.expr < b*2 | num.expr == 0)
#[1] 4911

no.expr.qq <- std.expr < b*2 | num.expr == 0
E <- E[!no.expr.qq,]

# Save nice gene-level text expression matrix for use with other tools
length(unique(featureData(E)$jetset.genesym)) == length(featureData(E)$jetset.genesym) # TRUE
M <- exprs(E)
rownames(M) <- featureData(E)$jetset.genesym
colnames(M) <- sub(".CEL.gz", "", colnames(M))

save(E,M,b, file="../E7307.may14.genelvl.exprfilt.RData")
write.table(M, file="../E7307.may14.genelvl.exprfilt.nice.tab", sep="\t", col.names=NA, row.names=TRUE, quote=FALSE)
