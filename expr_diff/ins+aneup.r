library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidygraph)
io = import('io')
sys = import('sys')
gset = import('data/genesets')
vp = import('../expr_cor/viper')
idmap = import('process/idmap')
util = import('./util')

args = sys$cmd$parse(
    opt('e', 'eset', 'gene expression RData', 'eset_Mad2PB.RData'),
    opt('f', 'config', 'yaml', '../config.yaml'),
    opt('c', 'cis', 'cis RData', '../cis_analysis/poisson.RData'),
    opt('i', 'ins', 'gene name of insert', 'Erg'),
    opt('n', 'network', 'aracne', '../data/networks/E-GEOD-13159.RData'),
    opt('o', 'outfile', 'results RData', 'ins/Erg.RData'),
    opt('p', 'plotfile', 'pdf', 'ins+aneup/Erg.pdf'),
    arg('sets', 'gene set .RData', arity='*',
        list.files("../data/genesets", "\\.RData", full.names=TRUE)))

om = function(x) idmap$orthologue(x, from="external_gene_name", to="mgi_symbol")
net = io$load(args$network) %>%
    mutate(Regulator = unname(om(Regulator)),
           Target = unname(om(Target))) %>%
    na.omit()
tf_net = filter(net, Target %in% Regulator)

cis = io$load(args$cis)$samples %>%
    filter(external_gene_name == args$ins) %>%
    select(sample, ins=external_gene_name)

dset = io$load(args$eset)
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

sets = io$load(args$sets) %>%
    setNames(tools::file_path_sans_ext(basename(args$sets))) %>%
    lapply(function(x) gset$filter(x, min=5, valid=rownames(eset)))

hl = io$read_yaml(args$config)$highlight_de

pdf(args$plotfile)
print(util$plot_pcs(idx, dset$pca, 1, 2, hl=cis$sample))

dviper = vp$diff_viper(dset$vs, net, eset$ins * eset$aneup0.2)
dcor = vp$diff_cor(dset$vs, tf_net, eset$ins * eset$aneup0.2)
print(vp$plot_subnet(dviper, dcor) + ggtitle("MI network"))

print(util$plot_volcano(res, hl) + ggtitle("gene"))
for (name in names(sets))
    print(util$plot_gset(res, sets[[name]]) + ggtitle(name))
dev.off()

save(res, file=args$outfile)
