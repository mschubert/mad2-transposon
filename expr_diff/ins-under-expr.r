library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidygraph)
io = import('io')
sys = import('sys')
gset = import('data/genesets')
vp = import('../link_cis_expr/cor_viper')
idmap = import('process/idmap')
rnaseq = import('process/rna-seq')
util = import('./util')

args = sys$cmd$parse(
    opt('e', 'eset', 'gene expression RData', 'eset_Mad2PB.RData'),
    opt('c', 'cis', 'cis RData', '../cis_analysis/poisson.RData'),
    opt('i', 'ins', 'gene name of insert', 'Erg'),
    opt('n', 'network', 'aracne', '../data/networks/E-GEOD-13159.RData'),
    opt('o', 'outfile', 'results RData', 'ins/Erg.RData'),
    opt('p', 'plotfile', 'pdf', 'ins/Erg.pdf'),
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
mgi = rownames(dset$vs) %>%
    idmap$gene(to="external_gene_name", dset="mmusculus_gene_ensembl")
idx = colData(eset) %>%
    as.data.frame() %>%
    left_join(cis) %>%
    mutate(expr = dset$vs[which(args$ins == mgi),],
           ins = ifelse(is.na(ins), 0, 1))
eset@colData = DataFrame(idx)
eset = eset[-which(args$ins == mgi),]

design(eset) = ~ tissue + type + ins + expr
res = DESeq2::estimateDispersions(eset) %>%
    DESeq2::nbinomLRT(reduced=~ tissue + type + ins, maxit=1000) %>%
    DESeq2::results() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("ensembl_gene_id") %>%
    mutate(gene_name = idmap$gene(ensembl_gene_id, from="ensembl_gene_id",
        to="external_gene_name", dset="mmusculus_gene_ensembl"))

sets = io$load(args$sets) %>%
    setNames(tools::file_path_sans_ext(basename(args$sets))) %>%
    lapply(function(x) gset$filter(x, min=5, valid=na.omit(mgi)))

pdf(args$plotfile)
print(util$plot_pcs(idx, dset$pca, 1, 2, hl=cis$sample))

#expr = assay(eset)
#rownames(expr) = idmap$gene(rownames(expr),
#    to="external_gene_name", dset="mmusculus_gene_ensembl")
dviper = vp$diff_viper(dset$vs, net, eset$ins)
dcor = vp$diff_cor(dset$vs, tf_net, eset$ins)
print(vp$plot_subnet(dviper, dcor) + ggtitle("MI network"))

print(util$plot_volcano(res) + ggtitle("gene"))
for (name in names(sets))
    print(util$plot_gset(res, sets[[name]]) + ggtitle(name))
dev.off()

save(res, file=args$outfile)
