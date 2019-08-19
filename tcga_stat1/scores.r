library(dplyr)
tcga = import('data/tcga')
idmap = import('process/idmap')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'sets.rds'),
    opt('o', 'outfile', 'rds', 'scores.rds'))

gsva_tcga = function(cohort, genes, sets) {
    expr = tcga$rna_seq(cohort, trans="raw")
    keep = rowMeans(expr, na.rm=TRUE) >= 10
    expr = tcga$rna_seq(cohort, trans="vst")
    rownames(expr) = idmap$gene(rownames(expr), to="external_gene_name")
    expr = expr[keep | rownames(expr) %in% genes,]
    expr = expr[!is.na(rownames(expr)) & rownames(expr) != "",]

    genes = expr[genes,,drop=FALSE]
    expr = expr[rownames(expr) != "STAT1",]
    gsva = GSVA::gsva(expr, sets, parallel.sz=1)
    scores = rbind(genes, gsva)
}

sets = readRDS(args$infile)
genes = c("STAT1", "PIAS1", "IFNG", "IL1B", "IFITM1",
          "CGAS", "TBK1", "IRF2", "IRF3", "TP53")
# don't have: "IFNA", "IFNB", "OAS1", "ISG54"

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

immune = readxl::read_xlsx("1-s2.0-S1074761318301213-mmc2.xlsx") %>%
    transmute(sample = `TCGA Participant Barcode`,
              aneup = as.numeric(`Aneuploidy Score`),
              stroma = as.numeric(`Stromal Fraction`),
              Lymphocytes = as.numeric(Lymphocytes),
              Lymph_Infiltration = as.numeric(`Lymphocyte Infiltration Signature Score`),
              Leuko_Fraction = as.numeric(`Leukocyte Fraction`),
              NK_activated = as.numeric(`NK Cells Activated`),
              NK_total = NK_activated + as.numeric(`NK Cells Resting`))

data = inner_join(scores, immune, by="sample") %>% na.omit()
saveRDS(data, file=args$outfile)
