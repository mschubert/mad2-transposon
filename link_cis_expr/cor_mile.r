io = import('io')
sys = import('sys')
util = import('./cor_util')

args = sys$cmd$parse(
    opt('s', 'select', 'yaml', 'interesting_sets.yaml'),
    opt('e', 'expr', 'expr RData', '../data/arrayexpress/E-GEOD-13159.RData'),
    opt('p', 'plotfile', 'pdf', 'cor_mile.pdf'),
    arg('genesets', 'RData files', arity='*',
        list.files("gsva_mile", "\\.RData$", full.names=TRUE))
)

dset = io$load(args$expr)
types = Biobase::pData(dset)$FactorValue..LEUKEMIA.CLASS.
expr = Biobase::exprs(dset)[c("ENSG00000134954", "ENSG00000157554"),]
rownames(expr) = c("ETS1", "ERG")

select = io$read_yaml(args$select)$expr_sets
sets = io$load(args$genesets)
sets = lapply(names(select), function(s) sets[[s]][select[[s]],,drop=FALSE])

mat = narray::stack(c(sets, list(expr)), along=1)
tmat = narray::split(mat, along=2, subsets=types)

pm = util$pcor(t(mat))

pdf(args$plotfile, 20, 15)
util$plot_cor_matrix(t(mat), text_color=NULL)
util$plot_pcor_net(pm, node_size=4, edge_size=1)
util$plot_bootstrapped_pcor(t(mat), node_size=4)

for (i in seq_along(tmat)) {
    name = names(tmat)[i]
    util$plot_cor_matrix(t(tmat[[i]]), title=name, text_color=NULL)
    try(util$plot_bootstrapped_pcor(t(tmat[[i]]), fdr=0.3, node_size=4, title=name))
}
dev.off()
