library(dplyr)
io = import('io')

mad2pb = io$load("../expr_diff/eset_Mad2PB.RData")
annot = as.data.frame(SummarizedExperiment::colData(mad2pb$eset))

mile = io$load("../expr_diff/eset_MILE.RData")#$expr # is there a coarse def?
meta = mile$meta

genes = c("Ets1", "Erg", "Stat1", "Pias1", "Stat3", "Pten", "Notch1")
c1 = cor(t(mad2pb$vs[genes, annot$type=="Myeloid"]))
c2 = cor(t(mad2pb$vs[genes, annot$type=="Other"]))
c3 = cor(t(mad2pb$vs[genes, annot$type=="T-cell"]))

genes = c("ETS1", "ERG", "STAT1", "PIAS1", "STAT3", "PTEN", "NOTCH1")
m1 = cor(t(mile$expr[genes, !is.na(meta$type) & meta$type=="Myeloid"]))
m2 = cor(t(mile$expr[genes, !is.na(meta$type) & meta$type=="B_like"]))
m3 = cor(t(mile$expr[genes, !is.na(meta$type) & meta$type=="T_ALL"]))

pdf("FigSx-corEtsErg.pdf", 6, 5)
corrplot::corrplot(c1, tl.cex=2, tl.col="black")
corrplot::corrplot(c2, tl.cex=2, tl.col="black")
corrplot::corrplot(c3, tl.cex=2, tl.col="black")

corrplot::corrplot(m1, tl.cex=2, tl.col="black")
corrplot::corrplot(m2, tl.cex=2, tl.col="black")
corrplot::corrplot(m3, tl.cex=2, tl.col="black")
dev.off()
