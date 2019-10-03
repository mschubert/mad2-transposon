library(dplyr)
library(ggplot2)
sys = import('sys')

args = sys$cmd$parse(
    opt('y', 'yaml', 'yaml', 'BT549_Stat1.yaml'),
    opt('e', 'expr', 'tsv', 'count_matrix_known_barcodes_STL_and_USS_genes.txt'),
    opt('o', 'outfile', 'rds', 'dset.rds'),
    opt('p', 'plotfile', 'pdf', 'dset.pdf'))

# map from yaml identifiers to count table identifiers
experiment = yaml::read_yaml(args$yaml)
samples = tibble(id=names(experiment$samples)) %>%
    mutate(genotype = sub("^([^_]+)_.*", "\\1", id),
           time = sub(".*_([0-9]+)[_-].*", "\\1", id),
           treatment = sub(".*(dmso|rev|ifng|0).*", "\\1", id),
           replicate = sapply(strsplit(id, "-"), function(x) x[[2]]))

# DESeq2 data set and size factors
dset = readr::read_tsv(args$expr)
counts = data.matrix(dset[-1])
rownames(counts) = sub("\\.[0-9]+", "", dset$gene_id)
colnames(counts) = samples$id
counts = counts[rowSums(counts) >= 10,]

# DESeq2 PCA of samples
eset = DESeq2::DESeqDataSetFromMatrix(counts, samples, ~1)
vst = DESeq2::varianceStabilizingTransformation(eset)

pdf(args$plotfile, 14, 8)
DESeq2::plotPCA(vst, intgroup=c("genotype", "time")) +
    ggrepel::geom_text_repel(aes(label=name))

DESeq2::plotPCA(vst[,!grepl("cgas_[^-]+-2", colnames(vst))], intgroup=c("genotype", "time")) +
    ggrepel::geom_text_repel(aes(label=name))
dev.off()

# save rds with count object
saveRDS(list(eset=eset, vst=vst, index=samples), file=args$outfile)
