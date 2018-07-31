library(BioNet)
library(DLBCL)
data(interactome)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('g', 'gene', 'gene-level poisson', '../cis_analysis/poisson.RData'),
    opt('d', 'diff_expr', 'aneuploidy DE', '../cis_analysis/aneup_de.RData'),
    opt('o', 'outfile', 'save results .RData', 'poisson_{set}.RData'),
    opt('p', 'plotfile', 'pdf', 'bionet.pdf'))

dset = io$load(args$gene)
samples = dset$samples
genes = dset$genes

stats = io$load(args$diff_expr) %>%
    mutate(score = -log10(p.value))

nodes = nodes(interactome) %>% sub("\\([0-9]+\\)", "", .)
scores = stats$score[match(nodes, toupper(stats$expr))]
scores[is.na(scores)] = 0
names(scores) = nodes(interactome)

diff_expr = stats$statistic[match(nodes, toupper(stats$expr))]
diff_expr[is.na(diff_expr)] = 0
names(diff_expr) = nodes(interactome)

module = runFastHeinz(interactome, scores)
plotModule(module, scores=scores, diff.expr=diff_expr)
dev.off()
