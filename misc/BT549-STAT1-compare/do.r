library(dplyr)
deseq = import('process/deseq')

# load nextflex, 2023 and 2024 data
# compare qPCR gene counts in seq data

# compare DMSO vs. 500 nM Rev 48 hours
c2021 = readr::read_tsv("../../data/rnaseq_stat1/count_matrix_known_barcodes_STL_and_USS_genes.txt")
counts = data.matrix(c2021[-1])
rownames(counts) = sub("\\.[0-9]+", "", c2021$gene_id)
experiment = yaml::read_yaml("../../data/rnaseq_stat1/BT549_Stat1.yaml")
samples = tibble(id=names(experiment$samples)) |>
    mutate(genotype = relevel(factor(sub("^([^_]+)_.*", "\\1", id)), "wt"),
           time = relevel(factor(sub(".*_([0-9]+)[_-].*", "\\1", id)), "0"),
           treatment = relevel(factor(sub(".*(dmso|rev|ifng|0).*", "\\1", id)), "dmso"),
           replicate = factor(sapply(strsplit(id, "-"), function(x) x[[2]])))
e2021 = DESeq2::DESeqDataSetFromMatrix(counts, samples, ~1)
e2021 = e2021[,grepl("wt_48_(rev|dmso)", e2021$id)]
r2021 = deseq$genes(e2021, ~ treatment) |>
    mutate(term = sub("^treatment", "r2021", term))

e2023 = readRDS("/DATA/m.schubert/nfcore-results/2023-03_smart3-BT549/rnaseq/star_salmon/salmon.merged.gene_counts.rds")
m2023 = readr::read_tsv("../../data/smart3seq/FF230302.tsv")
colData(e2023) = DataFrame(as.data.frame(m2023[match(e2023$names, m2023$sample),]))
e2023 = e2023[,e2023$genotype == "WT" & e2023$conc %in% c(0, 500)]
e2023$treatment = factor(e2023$treatment, levels=c("DMSO", "rev"))
rownames(colData(e2023)) = e2023$sample
assay(e2023) = data.matrix(assay(e2023))
r2023 = deseq$genes(DESeqDataSet(e2023, ~treatment)) |>
    mutate(term = sub("^treatment", "r2023", term))

e2024 = readRDS("/DATA/m.schubert/nfcore-results/2024-01_smart3-BT549/rnaseq/star_salmon/salmon.merged.gene_counts.rds")
#m2024 = readr::read_tsv("../../data/smart3seq/FF231221.tsv")
e2024 = e2024[,grepl("BT549_wt_(0|500)", e2024$names)]
e2024$cond = factor(sub("BT549_wt_([0-9]+_[0-9]+h).*", "\\1", e2024$names)) |> relevel("0_72h")
assay(e2024) = round(data.matrix(assay(e2024)))
r2024 = deseq$genes(DESeqDataSet(e2024, ~cond)) |>
    mutate(term = sub("^cond", "r2024", term))

dset = r2021 |> bind_rows(r2023) |> bind_rows(r2024) |>
    tidyr::unnest(genes) |>
    select(term, ensembl_gene_id, stat) |>
    tidyr::pivot_wider(names_from="term", values_from="stat")

# do corrplot all-vs-all
