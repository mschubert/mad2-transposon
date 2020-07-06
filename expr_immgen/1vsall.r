library(dplyr)
library(SummarizedExperiment)
library(ggplot2)
library(patchwork)
io = import('io')
plt = import('plot')

de_1vsall = function(eset, cur) {
    colData(eset)$is_cur = colData(eset)$Short == cur
    design(eset) = ~ is_cur
    res = DESeq2::DESeq(eset) %>%
        DESeq2::results(name="is_curTRUE") %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene_name") %>%
        as_tibble() %>%
        arrange(padj)
}

annot = readxl::read_xlsx("../data/immgen/immgen_keypop.xlsx")
immgen = readr::read_csv("../data/immgen/GSE109125_Gene_count_table.csv.gz")
coldata = tibble(id=colnames(immgen)[-1]) %>%
    mutate(Short = sub("#[0-9]+$", "", id)) #%>%
#    inner_join(annot)
expr = data.matrix(immgen[-1])
rownames(expr) = immgen$Gene_Symbol
expr = expr[,coldata$id]
expr = expr[rowMeans(expr) >= 1,]
# update gene symbols??
expr = DESeq2::DESeqDataSetFromMatrix(expr, coldata, ~1)
#immgen_de = lapply(unique(coldata$Short), de_1vsall, eset=expr)
immgen_de = clustermq::Q(de_1vsall, cur=unique(coldata$Short), const=list(eset=expr),
                         n_jobs=40, memory=2048, pkgs=c("dplyr", "SummarizedExperiment"))
names(immgen_de) = unique(coldata$Short)

tps = io$load("../expr_diff/eset_Mad2PB.RData")$eset
colData(tps)$Short = colnames(tps)
tps = tps[rowMeans(counts(tps)) >= 1,]
#tps_de = lapply(unique(colnames(tps)), de_1vsall, eset=tps)
tps_de = clustermq::Q(de_1vsall, cur=unique(colnames(tps)), const=list(eset=tps),
                      n_jobs=20, memory=2048, pkgs=c("dplyr", "SummarizedExperiment"))
names(tps_de) = unique(colnames(tps))

saveRDS(list(immgen_de=immgen_de, immgen_meta=annot,
             tps_de=tps_de, tps_meta=colData(tps)), file="1vsall.rds")

img = lapply(immgen_de, function(x) setNames(x$stat, x$gene_name)) %>%
    narray::stack(along=2)
img = img[rowMaxs(abs(img)) > 3 & !is.na(rownames(img)),]
tps = lapply(tps_de, function(x) setNames(x$stat, x$gene_name)) %>%
    narray::stack(along=2)
tps = tps[rowMaxs(abs(tps)) > 2 & !is.na(rownames(tps)),]
narray::intersect(img, tps, along=1) # loses a lot of genes

cosine = function(x, y) sum(x * y) / (sqrt(sum(x^2)) * sqrt(sum(y^2)))
res = narray::lambda(~ cosine(tps, img), along=c(tps=2, img=2))
clust = reshape2::melt(res) %>%
    plt$cluster(value ~ img + tps) %>%
    as_tibble()
pdf("1vsall.pdf", 20, 12)
ggplot(clust, aes(x=img, y=tps)) +
    geom_raster(aes(fill=value)) +
    scale_fill_distiller(palette="Spectral") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
dev.off()
