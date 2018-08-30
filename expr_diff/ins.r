library(DESeq2)
library(ggplot2)
library(dplyr)
io = import('io')
sys = import('sys')
de = import('./de')
set = import('./de_sets')

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'gene expression RData', 'de.RData'),
    opt('i', 'ins', 'gene name of insert', 'Erg'),
    opt('o', 'outfile', 'results RData', 'ins/Erg.RData'),
    opt('p', 'plotfile', 'pdf', 'ins/Erg.pdf'),
    arg('sets', 'gene set .RData', arity='*',
        list.files("../data/genesets", "\\.RData", full.names=TRUE)))

dset = io$load(args$diff_expr)
cis = dset$cis$samples %>%
    filter(external_gene_name == args$ins) %>%
    select(sample, ins=external_gene_name)
eset = dset$eset
idx = colData(eset) %>%
    as.data.frame() %>%
    left_join(cis) %>%
    mutate(ins = ifelse(is.na(ins), 0, 1))
eset@colData = DataFrame(idx)

design(eset) = ~ tissue + type + ins
res = DESeq2::estimateDispersions(eset) %>%
    DESeq2::nbinomWaldTest(maxit=1000) %>%
    de$extract_coef("ins") %>%
    mutate(gene_name = de$idmap$gene(ensembl_gene_id, from="ensembl_gene_id",
        to="external_gene_name", dset="mmusculus_gene_ensembl"))

sets = io$load(args$sets) %>%
    setNames(tools::file_path_sans_ext(basename(args$sets))) %>%
    lapply(function(x) set$gset$filter(x, min=5, max=200, valid=na.omit(res$gene_name)))

pdf(args$plotfile)
print(de$plot_pcs(idx, dset$pca, 1, 2, hl=cis$sample))
print(de$plot_volcano(res) + ggtitle("ins"))
for (name in names(sets))
    print(set$plot_gset(res, sets[[name]]) + ggtitle(name))
dev.off()

save(res, file=args$outfile)
