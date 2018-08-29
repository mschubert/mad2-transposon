sys = import('sys')
enr = import('tools/enrichr')
msdb = import('tools/msigdb')

args = sys$cmd$parse(
    opt('g', 'geneset', 'Identifier of the gene set'),
    opt('o', 'outfile', 'File to save gene set to'))

if (args$geneset %in% enr$dbs()$name)
    sets = enr$genes(args$geneset)
else if (args$geneset %in% msdb$dbs())
    sets = msdb$genes(args$geneset)
else
    stop("invalid gene set: ", args$geneset)

save(sets, file=args$outfile)
