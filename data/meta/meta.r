library(dplyr)
io = import('io')
idmap = import('process/idmap')

#TODO: RNA-derived tumor purity?
#TODO: what about B-cell genes?

rna = io$load("../rnaseq/assemble.RData")

meta = readr::read_tsv("meta.tsv")

genes = rna$counts %>%
    idmap$gene(to="external_gene_name", dset="mmusculus_gene_ensembl", summarize=sum)

tcr = genes[grepl("Trbv", rownames(genes)),] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::gather("id", "counts", -gene)
bcr = genes[grepl("Brbv", rownames(genes)),] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::gather("id", "counts", -gene)
icr = bind_rows(list(TCR=tcr, BCR=bcr), .id="type")

expression = genes[rownames(genes) %in% c("Mad2l1", "Trp53", "Msh2", "Pten"),] %>%
    edgeR::cpm(log=TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::gather("id", "vst", -gene)

weights = rna$idx %>%
    transmute(id=id, spleen=spleen_weight/1000, thymus=thymus_weight/1000) %>%
    tidyr::gather("tissue", "grams", -id)

save(meta, weights, icr, expression, file="meta.RData")
