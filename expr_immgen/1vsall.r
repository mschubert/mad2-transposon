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
#immgen = readr::read_csv("../data/immgen/GSE109125_Gene_count_table.csv.gz")
immgen = readr::read_csv("../data/immgen/GSE109125_Gene_count_table_GENCODE_vM25.csv.gz")
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
meta_tps = colData(tps) %>% as.data.frame() %>% as_tibble() %>%
    mutate(HighErg = counts(tps, normalized=TRUE)["Erg",] > 300)
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
meta = meta_tps %>%
    transmute(x1 = "type", x2="aneup", x3="high_erg",
              high_erg = HighErg,
              aneup = aneuploidy,
              tps = factor(sample, levels=levels(clust$tps)),
              type = type)
p1 = ggplot(meta, aes(x=x1, y=forcats::fct_reorder(tps, aneup))) +
    geom_raster(aes(fill=type)) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
p2 = ggplot(meta, aes(x=x2, y=forcats::fct_reorder(tps, aneup))) +
    geom_raster(aes(fill=aneup)) +
    scale_fill_distiller(palette="OrRd", direction=1) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
p3 = ggplot(meta, aes(x=x3, y=forcats::fct_reorder(tps, aneup))) +
    geom_raster(aes(fill=high_erg)) +
    scale_fill_manual(values=setNames(c("black", "white"), c(TRUE, FALSE))) +
    theme(axis.text.x = element_text(angle=45, hjust=1))
p4 = ggplot(clust %>% left_join(meta), aes(x=img, y=forcats::fct_reorder(tps, aneup))) +
    geom_raster(aes(fill=value)) +
    scale_fill_distiller(palette="Spectral") +
    theme(axis.text.x = element_text(angle=45, hjust=1))
pdf("1vsall.pdf", 20, 12)
p1 + p2 + p3 + p4 + plot_layout(widths=c(1,1,1,50), guides="collect")
dev.off()
