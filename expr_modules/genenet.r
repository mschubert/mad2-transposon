library(dplyr)
io = import('io')
seq = import('seq')
sys = import('sys')
gnet = import('tools/genenet')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression RData', '../data/rnaseq/assemble.RData'),
#    opt('c', 'copies', 'gene copy matrix', '../ploidy_compare/gene_copies.RData'),
    opt('c', 'cis', 'sites per sample', '../cis_analysis/poisson.RData'),
    opt('a', 'aneup', 'aneuploidy score', '../ploidy_compare/analysis_set.RData'),
    opt('o', 'outfile', 'results RData', 'aneup_de.RData'),
    opt('p', 'plotfile', 'pdf', 'aneup_de.pdf'))

aneup = io$load(args$aneup)

# adding this to network means 13 fewer samples, no assocs fdr<.25 w/ aneup
ins = io$load(args$cis)
is_cis = ins$result %>%
    filter(adj.p < 1e-6 | external_gene_name %in%
           c("Xrcc6", "Rapgef6", "Bcl11a", "Eed", "Suz12", "Arid2")) %>%
    pull(external_gene_name) %>%
    unique()
cis = ins$samples %>%
    filter(external_gene_name %in% is_cis) %>%
    narray::construct(n_ins ~ sample + external_gene_name, fill=0)
cis = t((!is.na(cis) & cis!=0) + 0)
rownames(cis) = paste("cis", rownames(cis), sep="_")

keep_gene = seq$coords$gene(dset="mmusculus_gene_ensembl", chromosomes=c(1:19,'X')) %>%
    pull(external_gene_name)

eset = io$load(args$expr) # filter on read count, variance?
rownames(eset$expr) = eset$gene
expr = eset$expr[eset$gene %in% keep_gene &
                 (rowMeans(eset$counts >= 10) >= 0.1 | eset$gene %in% cis),]
idx = eset$idx %>% transmute(sample=paste0(hist_nr, tissue), tissue, type)
narray::intersect(expr, idx$sample, aneup$sample, along=2) # cis: -13 samples
mat = rbind(aneup=aneup$aneup, expr)

# add IS as nodes?

# infer genenet
pcor = gnet$pcor(mat, fdr=0.25) # Ets1, Erg, aneup connected over 0.21 w/o CIS in net [aneup nothing with]

g = igraph::graph_from_data_frame(pcor, directed=FALSE) %>%
    as_tbl_graph() %>%
    tidygraph::mutate(community = as.factor(group_infomap()))

# subset w/ aneup+CIS cor [or only CIS for non-aneup related]
# plot this subset

# community detection?

# could also only take CIS genes for network
