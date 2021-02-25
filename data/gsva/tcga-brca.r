sys = import('sys')
tcga = import('data/tcga')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('s', 'setfile', 'RData', '../genesets/human/MSigDB_Hallmark_2020.rds'),
    opt('t', 'threads', 'num', '12'),
    opt('o', 'outfile', 'rds', 'tcga-brca/MSigDB_Hallmark_2020.rds')
)

sets = readRDS(args$setfile)
expr = tcga$rna_seq("BRCA", trans="vst")
rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")

scores = GSVA::gsva(expr, sets, parallel.sz=as.integer(args$threads))

saveRDS(scores, file=args$outfile)
