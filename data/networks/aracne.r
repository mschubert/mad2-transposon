library(dplyr)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')
gset = import('data/genesets')
ar = import('tools/aracne')

args = sys$cmd$parse(
    opt('i', 'infile', 'expression data', '../arrayexpress/E-GEOD-13159.RData'),
    opt('b', 'bootstraps', 'number', '100'),
    opt('o', 'outfile', 'save network to', 'E-GEOD-13159.RData'))

tfs = gset$go() %>%
    filter(id == "GO:0003700") %>%
    pull(hgnc_symbol)

expr = Biobase::exprs(io$load(args$infile))
rownames(expr) = unname(idmap$gene(rownames(expr), to="hgnc_symbol"))
expr = expr[!is.na(rownames(expr)),]

bs = as.integer(args$bootstraps)
clustermq::register_dopar_cmq(n_jobs=bs, memory=10240)
result = ar$aracne(expr, tfs, folder=".temp", bootstrap=bs)

save(result, file=args$outfile)
