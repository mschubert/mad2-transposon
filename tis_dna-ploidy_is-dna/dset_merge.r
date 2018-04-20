library(dplyr)
b = import('base')
io = import('io')

# load and merge aneuploidy scores
aneup = io$read_table("../ploidy_compare/compare_ploidy.tsv", header=TRUE) %>%
    select(-coverage, -tissue) %>%
    filter(!duplicated(data.frame(type, sample))) %>%
    tidyr::spread("type", "aneuploidy") %>%
    group_by(sample) %>%
    summarize(aneup = ifelse(sample == "",
        `WGS (single-cell)`,
        `WGS (merged)` %or% `WGS (30-cell)`)) %>%
    arrange(aneup) %>%
    filter(!sample %in% c("S", "T")) %>%
    mutate(sample = paste0(sub("[ST]", "", sample), tolower(sub("[0-9]+", "", sample))))

cis = io$load("../data/cis/cis_per_tumor.RData") %>%
    mutate(sample = ifelse(sample == "410bm", "410t", sample),
           sample = ifelse(sample == "415bm", "415s", sample),
           sample = ifelse(sample == "452s", "452t", sample)) %>%
    inner_join(aneup)

hits = cis %>%
    select(sample, gene_name, known_cancer, aneup, reads) %>%
    group_by(sample, gene_name, known_cancer, aneup) %>%
    summarize(reads = sum(reads)) %>%
    ungroup()

near = cis %>%
    select(sample, flanking, aneup, reads) %>%
    tidyr::unnest() %>%
    bind_rows(hits) %>%
    select(-ensembl_gene_id) %>%
    group_by(sample, gene_name, known_cancer, aneup) %>%
    summarize(reads = sum(reads)) %>%
    ungroup()

cancer_genes = near %>%
    filter(known_cancer) %>%
    pull(unique(gene_name))

cis_hits = narray::construct(reads ~ gene_name + sample, data=hits, fill=0)
cis_near = narray::construct(reads ~ gene_name + sample, data=near, fill=0)
cancer_hits = cis_hits[rownames(cis_hits) %in% cancer_genes,]
cancer_near = cis_near[rownames(cis_near) %in% cancer_genes,]

narray::intersect(aneup$sample, cis_hits, cis_near, cancer_hits, cancer_near, along=2)
save(aneup, cis_hits, cis_near, cancer_hits, cancer_near, file="dset.RData")
