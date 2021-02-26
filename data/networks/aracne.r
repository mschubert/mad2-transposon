library(dplyr)
sys = import('sys')
idmap = import('process/idmap')
gset = import('data/genesets')
ar = import('tools/aracne')

args = sys$cmd$parse(
    opt('i', 'infile', 'expression data', '../arrayexpress/E-GEOD-13159.rds'),
    opt('b', 'bootstraps', 'number', '100'),
    opt('o', 'outfile', 'save network to', 'aracne_E-GEOD-13159.rds')
)

tfs = gset$go(ontology="MF") %>%
    filter(go_id == "GO:0003700") %>% # DNA-binding transcription factor activity
    pull(hgnc_symbol)

expr = Biobase::exprs(readRDS(args$infile))
rownames(expr) = unname(idmap$gene(rownames(expr), to="hgnc_symbol"))
expr = expr[!is.na(rownames(expr)),]

bs = as.integer(args$bootstraps)
clustermq::register_dopar_cmq(n_jobs=bs, memory=10240)
result = ar$aracne(expr, tfs, folder=".temp", bootstrap=bs)

saveRDS(result, file=args$outfile)
