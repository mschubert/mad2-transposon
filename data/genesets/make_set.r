sys = import('sys')
idmap = import('process/idmap')
enr = import('tools/enrichr')
msdb = import('tools/msigdb')

args = sys$cmd$parse(
    opt('g', 'geneset', 'Identifier of the gene set', 'KEA_2015'),
    opt('o', 'outfile', 'File to save gene set to', '/dev/null'))

if (args$geneset %in% enr$dbs()$name)
    sets = enr$genes(args$geneset)
else if (args$geneset %in% msdb$dbs())
    sets = msdb$genes(args$geneset)
else
    stop("invalid gene set: ", args$geneset)

mouse = stack(sets)
mouse$values = unname(idmap$orthologue(mouse$values, from="external_gene_name", to="mgi_symbol"))
mouse = unstack(na.omit(mouse))

save(mouse, file=args$outfile)
