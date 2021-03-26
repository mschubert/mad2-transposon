library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidygraph)
sys = import('sys')
gset = import('genesets')
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
    opt('p', 'plotfile', 'pdf', 'ins+aneup/Erg.pdf'),
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
idx = colData(eset) %>%
    as.data.frame() %>%
    left_join(cis) %>%
    mutate(ins = ifelse(is.na(ins), 0, 1),
           aneup0.2 = pmin(aneuploidy, 0.2))
eset@colData = DataFrame(idx)

design(eset) = ~ tissue + type + type:aneup0.2 + ins + ins:aneup0.2
res = DESeq2::estimateDispersions(eset) %>%
    DESeq2::nbinomLRT(reduced=~ tissue + type + type:aneup0.2 + ins, maxit=1000) %>%
    DESeq2::results() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene_name") %>%
    arrange(padj, pvalue)

sets = lapply(args$sets, readRDS) %>%
    setNames(tools::file_path_sans_ext(basename(args$sets))) %>%
    lapply(function(x) gset$filter(x, min=5, valid=rownames(eset)))

hl = yaml::read_yaml(args$config)$highlight_de

pdf(args$plotfile)
print(util$plot_pcs(idx, dset$pca, 1, 2, hl=cis$sample))

dviper = vp$diff_viper(dset$vs, net, eset$ins * eset$aneup0.2)
dcor = vp$diff_cor(dset$vs, tf_net, eset$ins * eset$aneup0.2)
print(vp$plot_subnet(dviper, dcor) + ggtitle("MI network"))

print(util$plot_volcano(res, hl) + ggtitle("gene"))
for (name in names(sets))
    print(util$plot_gset(res, sets[[name]]) + ggtitle(name))
dev.off()

saveRDS(res, file=args$outfile)
