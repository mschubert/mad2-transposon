library(dplyr)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression', '../data/rnaseq/assemble.RData'),
    opt('s', 'geneset', 'gene set RData', '../data/genesets/GO_Biological_Process_2018.RData'),
    opt('o', 'outfile', 'save RData', 'gsva_mad2pb/GO_Biological_Process_2018.RData')
)

sets = io$load(args$geneset)

dset = io$load(args$expr)
if (is.list(dset)) {
    expr = dset$expr
    rownames(expr) = idmap$gene(rownames(expr), to="mgi_symbol", dset="mmusculus_gene_ensembl")
    expr = expr[!is.na(rownames(expr)),]
} else {
    expr = Biobase::exprs(dset)
    rownames(expr) = idmap$orthologue(rownames(expr), to="mgi_symbol")
    expr = expr[!is.na(rownames(expr)),]
}

scores = GSVA::gsva(expr, sets, parallel.sz=1)

save(scores, file=args$outfile)
