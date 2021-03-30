library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidygraph)
sys = import('sys')
gset = import('data/genesets')
vp = import('../expr_cor/viper')
idmap = import('process/idmap')
util = import('./util')

args = sys$cmd$parse(
    opt('e', 'eset', 'gene expression rds', 'eset_Mad2PB.rds'),
    opt('f', 'config', 'yaml', '../config.yaml'),
    opt('c', 'cis', 'cis rds', '../cis_analysis/poisson.rds'),
    opt('i', 'ins', 'gene name of insert', 'Erg'),
    opt('n', 'network', 'aracne', '../data/networks/E-GEOD-13159.rds'),
    opt('o', 'outfile', 'results rds', 'ins/Erg.rds'),
    opt('p', 'plotfile', 'pdf', 'ins/Erg.pdf'),
    arg('sets', 'gene set .rds', arity='*',
        list.files("../data/genesets", "\\.rds", full.names=TRUE))
)

om = function(x) idmap$orthologue(x, from="external_gene_name", to="mgi_symbol")
net = readRDS(args$network) %>%
    mutate(Regulator = unname(om(Regulator)),
           Target = unname(om(Target))) %>%
    na.omit()
tf_net = filter(net, Target %in% Regulator)

cis = readRDS(args$cis)$samples %>%
    filter(external_gene_name == args$ins) %>%
    select(sample, ins=external_gene_name)

dset = readRDS(args$eset)
eset = dset$eset
eset$ins = ifelse(eset$sample %in% cis$sample, 1, 0)
eset$aneup0.2 = pmin(eset$aneuploidy, 0.2)
idx = as.data.frame(colData(eset))

design(eset) = ~ tissue + type * aneup0.2 + ins
res = DESeq2::DESeq(eset) %>%
    DESeq2::results(name="ins") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_name") %>%
    as_tibble() %>%
    arrange(padj, pvalue)
#design(eset) = ~ tissue + type + aneup0.2 + ins
#res = DESeq2::estimateDispersions(eset) %>%
#    DESeq2::nbinomLRT(reduced=~ tissue + type + aneup0.2, maxit=1000) %>%
#    DESeq2::results() %>%
#    as.data.frame() %>%
#    tibble::rownames_to_column("gene_name") %>%
#    mutate(stat = sign(log2FoldChange) * stat)

sets = lapply(args$sets, readRDS) %>%
    setNames(tools::file_path_sans_ext(basename(args$sets))) %>%
    lapply(function(x) gset$filter(x, min=5, valid=rownames(eset)))

hl = yaml::read_yaml(args$config)$highlight_de

pdf(args$plotfile)
print(util$plot_pcs(idx, dset$pca, 1, 2, hl=cis$sample))

dviper = vp$diff_viper(dset$vs, net, eset$ins)
dcor = vp$diff_cor(dset$vs, tf_net, eset$ins)
print(vp$plot_subnet(dviper, dcor) + ggtitle("MI network"))

print(util$plot_volcano(res, hl) + ggtitle("gene"))
for (name in names(sets))
    print(util$plot_gset(res, sets[[name]]) + ggtitle(name))
dev.off()

saveRDS(res, file=args$outfile)
