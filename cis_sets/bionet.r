library(BioNet)
library(DLBCL)
library(ggraph)
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
stats = dset$result

nodes = nodes(interactome) %>% sub("\\([0-9]+\\)", "", .)
scores = stats$p.value[match(nodes, toupper(stats$external_gene_name))]
scores = -log10(scores)
scores[is.na(scores)] = 0
scores = pmax(scores - 2, 0) # only keep p<1e-x
names(scores) = nodes(interactome)

#bum = fitBumModel(scores, plot=TRUE)
#scores = scoreNodes(network=interactome, fb=bum, fdr=1e-10)
module = runFastHeinz(interactome, scores) %>%
    igraph::igraph.from.graphNEL()

ggraph(module) +
    geom_edge_link(alpha=0.2) +
    geom_node_point(alpha=0.2) +
    geom_node_text(aes(label = sub("^([A-Z0-9]+).*", "\\1", name)), size=2, repel=TRUE)

dev.off()
