library(dplyr)
io = import('io')
idmap = import('process/idmap')

#TODO: RNA-derived tumor purity?
#TODO: what about B-cell genes?

rna = io$load("../rnaseq/assemble.RData")

genes = edgeR::cpm(rna$counts, log=FALSE) %>%
    idmap$gene(to="external_gene_name", dset="mmusculus_gene_ensembl", summarize=sum)
genesT = genes[grepl("Trbv", rownames(genes)),] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::gather("id", "cpm", -gene)
genesB = genes[grepl("Brbv", rownames(genes)),] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::gather("id", "cpm", -gene)
other = genes[rownames(genes) %in% c("Mad2l1", "Trp53", "Msh2", "Pten"),] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::gather("id", "cpm", -gene)

meta = rna$idx %>% select(id, type=`Tumour type`, stage=`Early/ Late`)
gene = bind_rows(list(TCR=genesT, BCR=genesB, General=other), .id="type") %>%
    select(id, everything())
weights = rna$idx %>%
    select(id, spleen=spleen_weight, thymus=thymus_weight) %>%
    tidyr::gather("tissue", "mg", -id)

save(meta, weights, gene, file="meta.RData")
