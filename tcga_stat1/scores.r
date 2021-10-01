library(dplyr)
tcga = import('data/tcga')
idmap = import('process/idmap')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'sets/stat1TGs.rds'),
    opt('o', 'outfile', 'rds', 'scores/stat1TGs.rds')
)

gsva_tcga = function(cohort, genes, sets) {
    message(cohort)
    expr = tcga$rna_seq(cohort, trans="raw")
    keep = rowMeans(expr, na.rm=TRUE) >= 10
    expr = tcga$rna_seq(cohort, trans="vst")
    rownames(expr) = idmap$gene(rownames(expr), to="external_gene_name")
    expr = expr[keep | rownames(expr) %in% genes,]
    expr = expr[!is.na(rownames(expr)) & rownames(expr) != "", !duplicated(colnames(expr))]

    genes = expr[intersect(rownames(expr), genes),,drop=FALSE]
    if (length(sets) == 0)
        return(genes)

    expr = expr[! rownames(expr) %in% genes,]
    gsva = GSVA::gsva(expr, sets, parallel.sz=1)
    scores = rbind(genes, gsva)
}

dset = readRDS(args$infile)
sets = dset$sets
genes = dset$genes

cohorts = c("BRCA", "LUAD", "COAD", "PRAD", "SKCM")
gsva = lapply(cohorts, gsva_tcga, sets=sets, genes=genes)
scores = gsva %>%
    narray::stack(along=2) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    mutate(cohort = tcga$barcode2study(sample)) %>%
    filter(substr(sample, 14, 16) == "01A") %>%
    mutate(sample = substr(sample, 1,12))
#    tcga$filter(cancer==TRUE, primary==TRUE)

aneup = lapply(cohorts, tcga$aneuploidy) %>%
    bind_rows() %>%
    inner_join(tcga$purity() %>% select(Sample, purity=estimate)) %>%
    transmute(sample=substr(Sample,1,12), aneup=aneup_log2seg / purity)

immune = readxl::read_xlsx("1-s2.0-S1074761318301213-mmc2.xlsx") %>%
    transmute(sample = `TCGA Participant Barcode`,
              stroma = as.numeric(`Stromal Fraction`),
              Lymphocytes = as.numeric(Lymphocytes),
              Lymph_Infiltration = as.numeric(`Lymphocyte Infiltration Signature Score`),
              Leuko_Fraction = as.numeric(`Leukocyte Fraction`),
              Macrophages = as.numeric(`Macrophage Regulation`),
              NK_activated = as.numeric(`NK Cells Activated`),
              NK_total = NK_activated + as.numeric(`NK Cells Resting`)) %>%
    inner_join(aneup)

data = inner_join(scores, immune, by="sample") %>% na.omit()
saveRDS(data, file=args$outfile)
