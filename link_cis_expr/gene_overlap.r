io = import('io')
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
    opt('s', 'select', 'yaml which sets to use', './interesting_sets.yaml'),
    opt('p', 'plotfile', 'pdf', 'gene_overlap.pdf'),
    arg('genesets', 'RData files', arity='*',
        list.files("../data/genesets", "\\.RData$", full.names=TRUE))
)

select = io$read_yaml(args$select)
sets = io$load(args$genesets)

pdf(args$plotfile)
dev.off()
