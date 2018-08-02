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
cis = io$load(args$cis)
is_cis = cis$result %>% filter(adj.p < 1e-3) %>% pull(external_gene_name) %>% unique()

eset = io$load(args$expr) # filter on read count, variance?
rownames(eset$expr) = eset$gene
expr = eset$expr[rowMeans(eset$counts >= 10) >= 0.1 | eset$gene %in% cis,]
idx = eset$idx %>% transmute(sample=paste0(hist_nr, tissue), tissue, type)
narray::intersect(expr, idx$sample, aneup$sample, along=2)
mat = rbind(aneup=aneup$aneup, expr)

# add IS as nodes?

# infer genenet
pcor = gnet$pcor(mat, fdr=0.1)

# subset w/ aneup+CIS cor [or only CIS for non-aneup related]
# plot this subset

# community detection?

# could also only take CIS genes for network
