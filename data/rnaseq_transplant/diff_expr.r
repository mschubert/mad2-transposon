library(dplyr)
library(DESeq2)
sys = import('sys')
idmap = import('process/idmap')
gset = import('genesets')

de_comparison = function(rec, sets) {
    keep = as.data.frame(rec$samples) %>% mutate(in_comparison=1) %>%
        right_join(as.data.frame(colData(eset))) %>%
        filter(!is.na(in_comparison))

    eset2 = eset[,keep$Sample_ID]
    design(eset2) = as.formula(rec$design)

    genes = DESeq2::DESeq(eset2) %>%
        DESeq2::results(name=rec$extract) %>%
        as.data.frame() %>% tibble::rownames_to_column("ensembl_gene_id") %>%
        as_tibble() %>%
        mutate(label = idmap$gene(ensembl_gene_id, to="external_gene_name")) %>%
        select(ensembl_gene_id, label, everything()) %>%
        arrange(padj, pvalue)

    sres = lapply(sets, function(s) gset$test_lm(genes, s, add_means="log2FoldChange"))
    c(list(genes=genes), sres)
}

args = sys$cmd$parse(
    opt('e', 'eset', 'rds', 'samples.rds'),
    opt('c', 'comps', 'yaml', 'de_clines.yaml'),
    opt('o', 'outfile', 'rds', 'de_clines.rds'),
    opt('x', 'xlsfile', 'xlsx', 'de_clines.xlsx')
)

comps = yaml::read_yaml(args$comps)$comparisons
eset = readRDS(args$eset)

sets = gset$get_mouse(c(
    "MSigDB_Hallmark_2020",
    "GO_Biological_Process_2021",
    "DoRothEA",
    "CIN"
), conf="A")

res = lapply(comps, de_comparison, sets=sets)
saveRDS(res, file=args$outfile)

genes = lapply(res, function(x) x$genes)
writexl::write_xlsx(genes, args$xlsfile)
