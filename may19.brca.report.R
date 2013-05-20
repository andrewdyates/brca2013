load("../BRCA.R.may17.RData")
load("~/tftargets/tf.adj.RData")
load("~/pina/pina.adj.RData")
load("../BRCA.CLS.RData")
load("../BRCA1.dcor.RData")
load("../BRCA2.dcor.RData")
GLYPH.COLS <- c("#ffffff", "#40a185", "#2688bf", "#5b51a5", "#000000", "#a00d42", "#d7424c", "#eb6532", "#ffffff")

# Get PCC
load("../E7307.may14.genelvl.exprfilt.RData")
BRCA1.PCC <- apply(M, 1, function(x) cor(x,M[which(rownames(M) == "BRCA1"),]))
save(BRCA1.PCC, file="../BRCA1.pcc.RData")
BRCA2.PCC <- apply(M, 1, function(x) cor(x,M[which(rownames(M) == "BRCA2"),]))
save(BRCA2.PCC, file="../BRCA2.pcc.RData")

# 1: get transcription factors targeting BRCA1,2
colnames(TF.ADJ)[which(TF.ADJ[which(rownames(TF.ADJ)=="BRCA1"),]==1)]
colnames(TF.ADJ)[which(TF.ADJ[which(rownames(TF.ADJ)=="BRCA2"),]==1)]
brca2.tf <- colnames(TF.ADJ)[which(TF.ADJ[which(rownames(TF.ADJ)=="BRCA2"),]==1)]
split(BRCA2.CLS[names(BRCA2.CLS) %in% brca2.tf], BRCA2.CLS[names(BRCA2.CLS) %in% brca2.tf])

# 2: enrichment plots
plot.enrichments <- function(tab, PCC, bb, name) {
  # ppi
  ppis <- PINA.ADJ[which(rownames(PINA.ADJ)==name),]
  ppis.syms <- names(ppis)[ppis==1]
  pcc.ppi <- cumsum(names(sort(PCC,decreasing=T)) %in% ppis.syms)
  rand.ppi <- cumsum(sample(names(PCC)) %in%  ppis.syms)

  pdf(paste0("ppi.enrich.",name,".pdf"))
  plot(fit.v(cumsum(names(sort(bb,decreasing=T)) %in% ppis.syms)), col="#ff0000", type="l", ylab="PPI Hits")
  lines(fit.v(rand.ppi), type="l", col="#999999")
  lines(fit.v(cumsum(tab$"1"$ppi)), col=GLYPH.COLS[2])
  lines(fit.v(cumsum(tab$"2"$ppi)), col=GLYPH.COLS[3])
  lines(fit.v(cumsum(tab$"3"$ppi)), col=GLYPH.COLS[4])
  lines(fit.v(cumsum(tab$"4"$ppi)), col=GLYPH.COLS[5])
  lines(fit.v(cumsum(tab$"5"$ppi)), col=GLYPH.COLS[6])
  lines(fit.v(pcc.ppi), col="#ffff00")
  dev.off()

  # transcription factors
  tfs <- colnames(TF.ADJ)
  tfs.targ <- TF.ADJ[which(rownames(TF.ADJ)==name),]
  tfs.tags.syms <- names(tfs.targ)[tfs.targ==1]

  pcc.tf <- cumsum(names(sort(PCC[names(PCC) %in% tfs], decreasing=T)) %in% tfs.tags.syms)
  rand.tf <- cumsum(sample(names(PCC)[names(PCC) %in% tfs]) %in% tfs.tags.syms)
  pdf(paste0("tf.enrich.",name,".pdf"))
  plot(fit.v(pcc.tf), type="l", col="#ffff00", ylab="TF Targets Hits", xlab="Rank of TF genes only")
  lines(fit.v(rand.tf), type="l", col="#999999")
  lines(fit.v(cumsum(tab$"1"$trans.targ.v[rownames(tab$"1") %in% tfs])), col=GLYPH.COLS[2])
  lines(fit.v(cumsum(tab$"2"$trans.targ.v[rownames(tab$"2") %in% tfs])), col=GLYPH.COLS[3])
  lines(fit.v(cumsum(tab$"3"$trans.targ.v[rownames(tab$"3") %in% tfs])), col=GLYPH.COLS[4])
  lines(fit.v(cumsum(tab$"4"$trans.targ.v[rownames(tab$"4") %in% tfs])), col=GLYPH.COLS[5])
  lines(fit.v(cumsum(tab$"5"$trans.targ.v[rownames(tab$"5") %in% tfs])), col=GLYPH.COLS[6])
  dev.off()

  # PCC by class
  pdf(paste0("pcc.per.class.",name,".pdf"))
  plot(fit.v(sort(PCC,decreasing=T)), type="l", col="#ffff00", ylab="PCC")
  lines(fit.v(PCC[match(rownames(tab$"1"),names(PCC))]), col=GLYPH.COLS[2])
  lines(fit.v(PCC[match(rownames(tab$"2"),names(PCC))]), col=GLYPH.COLS[3])
  lines(fit.v(PCC[match(rownames(tab$"3"),names(PCC))]), col=GLYPH.COLS[4])
  lines(fit.v(PCC[match(rownames(tab$"4"),names(PCC))]), col=GLYPH.COLS[5])
  dev.off()
}
fit.v <- function(v, n=1000) {
  if (length(v) > n)
    v <- v[1:n]
  v
}
b1 <- (BRCA1.CLS==1 & BRCA1.DCOR>=0.7) | BRCA1.CLS %in% c(2,3) & BRCA1.DCOR >= 0.2
b2 <- (BRCA2.CLS==1 & BRCA2.DCOR>=0.7) | BRCA2.CLS %in% c(2,3) & BRCA2.DCOR >= 0.2
plot.enrichments(BRCA1.R$table, BRCA1.PCC, BRCA1.DCOR[b1], "BRCA1")
plot.enrichments(BRCA2.R$table, BRCA2.PCC, BRCA2.DCOR[b2], "BRCA2")



# 3: fraction per class for PPI, TF, TFT, ALL

