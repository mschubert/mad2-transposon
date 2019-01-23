library(dplyr)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression', '../data/arrayexpress/E-GEOD-13159.RData'),
    opt('s', 'geneset', 'gene set RData', '../data/genesets/human/GO_Biological_Process_2018.RData'),
    opt('o', 'outfile', 'save RData', 'gsva_mad2pb/GO_Biological_Process_2018.RData')
)

sets = io$load(args$geneset)
dset = io$load(args$expr)
expr = Biobase::exprs(dset)
rownames(expr) = idmap$gene(rownames(expr), to="external_gene_name")
expr = expr[!is.na(rownames(expr)),]

scores = GSVA::gsva(expr, sets, parallel.sz=1)
save(scores, file=args$outfile)
