library(dplyr)
library(DESeq2)
library(ggplot2)
library(patchwork)
sys = import('sys')
gset = import('genesets')
idmap = import('process/idmap')
plt = import('plot')

plot_pca = function(vst) {
    pcadata = DESeq2::plotPCA(vst, intgroup=c("Sample_ID", "Sample_name",
        "driver", "STAT1", "KIF2C/dnMCAK", "Sample_type"), returnData=TRUE)
    pcavar = round(100 * attr(pcadata, "percentVar"))

    ggplot(pcadata, aes(x=PC1, y=PC2)) +
        geom_point(aes(color=driver, shape=STAT1, size=KIF2C.dnMCAK), alpha=0.9) +
        ggrepel::geom_text_repel(aes(label=Sample_name)) +
        scale_size_manual(values=c(KIF2C=5, dnMCAK=10)) +
        xlab(paste0("PC1: ", pcavar[1], "% variance")) +
        ylab(paste0("PC2: ", pcavar[2], "% variance"))
}

plot_gsva = function(df, set_name) {
    df %>% filter(set == set_name) %>%
        ggplot(aes(x=Sample_name, y=value)) +
            geom_point(aes(color=driver, shape=STAT1, size=`KIF2C/dnMCAK`), alpha=0.9) +
            theme(axis.text.x = element_text(angle=30, hjust=1)) +
            facet_wrap(~ Sample_type, scales="free_x") +
            ggtitle(set_name) +
            scale_size_manual(values=c(KIF2C=5, dnMCAK=10)) +
            scale_shape_manual(values=c(WT=16, KO=25))
}

plot_HMpca = function(scores) {
    meta2 = meta[match(meta$Sample_ID, colnames(scores)),]
    pc1 = prcomp(t(scores), scale.=TRUE)
    plt$pca(pc1, aes(x=PC1, y=PC2), meta2, biplot=TRUE, bi_size=2.5, bi_arrow=0.05) +
        geom_point(aes(color=driver, shape=STAT1, size=`KIF2C/dnMCAK`), alpha=0.9) +
        ggrepel::geom_text_repel(aes(label=Sample_name), size=3.5)
}

args = sys$cmd$parse(
    opt('i', 'infile', 'counts', 'count_matrix_known_barcodes_STL_and_USS_genes.txt.gz'),
    opt('s', 'samples', 'tsv', 'samples.tsv'),
    opt('o', 'outfile', 'rds', 'samples.rds'),
    opt('p', 'plotfile', 'pdf', 'samples.pdf')
)

meta = readr::read_tsv(args$samples) %>%
    mutate(driver = paste0(Myc, p53),
           driver = sub("OEXWT", "Myc", driver, fixed=TRUE),
           driver = sub("WTKO", "p53", driver, fixed=TRUE))
meta$Sample_ID = paste0("SU_", meta$Sample_ID)
reads = readr::read_tsv(args$infile)
stopifnot(meta$Sample_ID == colnames(reads)[-1])

eset = DESeq2::DESeqDataSetFromMatrix(reads[-1], meta, ~1)
rownames(eset) = sub("\\.[0-9]+$", "", reads$gene_id)
vs = DESeq2::varianceStabilizingTransformation(eset)
rownames(vs) = idmap$gene(rownames(eset), to="hgnc_symbol")

sets = gset$get_mouse("MSigDB_Hallmark_2020")
scores = GSVA::gsva(assay(vs), sets)
scores = rbind(scores, assay(vs)[c("Myc", "Trp53", "Stat1"),])
names(dimnames(scores)) = c("set", "Sample_ID")

df = inner_join(meta, reshape2::melt(scores))

pdf(args$plotfile, 12, 9)

plot_pca(vs) + ggtitle("PCA all samples")
plot_pca(vs[,eset$Sample_type == "Cell line"]) + ggtitle("PCA cell lines")
plot_pca(vs[,eset$Sample_type == "Tumor"]) + ggtitle("PCA tumors")

plot_gsva(df, "Trp53") /
plot_gsva(df, "p53 Pathway") /
plot_gsva(df, "IL-6/JAK/STAT3 Signaling") + plot_layout(guides="collect")

plot_gsva(df, "Stat1") /
plot_gsva(df, "Interferon Alpha Response") /
plot_gsva(df, "Interferon Gamma Response") + plot_layout(guides="collect")

plot_gsva(df, "Myc") /
plot_gsva(df, "Myc Targets V1") /
plot_gsva(df, "Myc Targets V2") + plot_layout(guides="collect")

plot_HMpca(scores) + ggtitle("Hallmark PCA all samples")
plot_HMpca(scores[,eset$Sample_type == "Cell line"]) + ggtitle("Hallmark PCA cell lines")
plot_HMpca(scores[,eset$Sample_type == "Tumor"]) + ggtitle("Hallmark PCA tumors")

dev.off()
