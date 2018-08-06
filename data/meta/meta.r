library(dplyr)
io = import('io')
idmap = import('process/idmap')

#TODO: RNA-derived tumor purity?

rna = io$load("../rnaseq/assemble.RData")

meta = readr::read_tsv("meta.tsv")

genes = idmap$gene(c("Mad2l1", "Trp53", "Msh2", "Pten"), from="external_gene_name",
                   to="ensembl_gene_id", dset="mmusculus_gene_ensembl")

cd = data.frame(id=colnames(rna$counts))
eset = DESeq2::DESeqDataSetFromMatrix(rna$counts, colData=cd, ~1) %>%
    DESeq2::estimateSizeFactors() %>%
    DESeq2::counts(normalized=TRUE)

eset = eset[genes,]
rownames(eset) = names(genes)
expression = as.data.frame(eset) %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::gather("id", "vst", -gene)

save(meta, expression, file="meta.RData")
