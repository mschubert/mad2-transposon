library(dplyr)
library(ggplot2)
theme_set(cowplot::theme_cowplot())
b = import('base')
st = import('stats')
idmap = import('process/idmap')
sys = import('sys')
util = import('./util')

args = sys$cmd$parse(
    opt('m', 'mile', 'expr rdata', '../data/arrayexpress/E-GEOD-13159.rds'),
    opt('a', 'aneup', 'mile aneup rdata', '../misc/subtype_overlay/aneup_scores_mile.rds'),
    opt('o', 'outfile', 'rds', 'eset_MILE.rds'),
    opt('p', 'plotfile', 'pdf', 'eset_MILE.pdf')
)

aneup = readRDS(args$aneup)$aneup
mile = readRDS(args$mile)
expr = Biobase::exprs(mile)

rownames(expr) = idmap$gene(rownames(expr), to="external_gene_name")
expr = expr[!is.na(rownames(expr)) & !duplicated(rownames(expr)),]

meta = Biobase::pData(mile) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("id") %>%
    select(id, annot=FactorValue..LEUKEMIA.CLASS.) %>%
    mutate(lineage = case_when(
               grepl("[AC]ML", annot) ~ "Myeloid",
               grepl("[AC]LL", annot) ~ "Lymphoid",
               TRUE ~ as.character(NA)
           ),
           type = case_when(
               lineage == "Myeloid" ~ "Myeloid",
#               annot == "mature B-ALL with t(8;14)" ~ "B_ALL",
               annot == "T-ALL" ~ "T_ALL",
               lineage == "Lymphoid" & annot != "CLL" ~ "B_like"
           ),
           aneuploidy = aneup$aneuploidy[match(.$id, aneup$sample)])

idx = meta %>% mutate(sample=NA, tissue=lineage)
narray::intersect(meta$id, idx$id, expr, along=2)
pca = prcomp(t(expr), center=TRUE, scale=FALSE)

pdf(args$plotfile)
print(util$plot_pcs(idx, pca, 1, 2))
print(util$plot_pcs(idx, pca, 3, 4))
print(util$plot_pcs(idx, pca, 5, 6))
dev.off()

saveRDS(list(expr=expr, meta=meta, pca=pca), file=args$outfile)
