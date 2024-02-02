library(dplyr)
library(DESeq2)
library(ggplot2)
gset = import('genesets')

dset = readr::read_tsv("zhenya_redo/Floris_240126.all_counts")
dmat = data.matrix(dset[-1])
rownames(dmat) = dset[[1]]
colnames(dmat) = c("DMSO72h", "Reversine500nM72h", "Stat1KO.Reversine500nM72h")

eset = DESeqDataSetFromMatrix(round(dmat), ~1, colData=data.frame(sample=smp)) |>
    estimateSizeFactors()
vs = DESeq2::varianceStabilizingTransformation(eset)

qPCR = counts(eset, normalized=TRUE)[c("CXCL10", "ISG15", "IFNB1", "IL6", "CXCL8", "CCL5"),] |>
    reshape2::melt()
ggplot(qPCR, aes(x=Var2, fill=Var2, y=value)) +
    facet_wrap(~Var1) +
    geom_col() +
    scale_y_continuous(trans="log1p", breaks=c(1,5,20,100,500)) +
    theme(axis.text.x = element_text(angle=90, hjust=1))

sets = gset$get_human("MSigDB_Hallmark_2020")
scores = GSVA::gsva(assay(vs), sets)
pdf("Rplots.pdf",5,10)
ComplexHeatmap::Heatmap(scores, row_dend_reorder=TRUE)
dev.off()
