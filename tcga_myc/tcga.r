library(dplyr)
library(SummarizedExperiment)
sys = import('sys')
idmap = import('process/idmap')
tcga = import('data/tcga')

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'tcga.rds')
)

eset = tcga$rna_seq("BRCA", trans="vst", annot=TRUE)
expr = assay(eset)
rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
expr = expr[rownames(expr) != "" & !is.na(rownames(expr)) & apply(expr, 1, sd) > 0.5,]
purity = tcga$purity() %>%
    transmute(Sample = Sample, purity=estimate)
meta = colData(eset) %>%
    as.data.frame() %>%
    as_tibble() %>%
    transmute(Sample = sample,
              tumor_stage = sub("stage ", "", tumor_stage),
#              tumor_grade = tumor_grade, # not reported
              OS_time = pmax(days_to_death, days_to_last_follow_up, na.rm=TRUE),
              vital_status = factor(vital_status, levels=c("alive", "dead")),
#              site_orig = tissue_or_organ_of_origin, # c22.0
              sex = gender,
              year_of_birth = year_of_birth) %>%
    left_join(purity) %>%
    tcga$filter(along="Sample", cancer=TRUE, primary=TRUE, vial="A")
cna = tcga$cna_genes("BRCA", gene="hgnc_symbol")
narray::intersect(expr, cna, along=1)
narray::intersect(expr, cna, meta$Sample, along=2)

mut = tcga$mutations() %>%
    filter(Study == "BRCA",
           Variant_Classification != "Silent") %>%
    select(-Study)

dset = list(meta=meta, expr=expr, mut=mut, cna=cna)
saveRDS(dset, file=args$outfile)
