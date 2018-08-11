library(BioNet)
library(DLBCL)
library(ggraph)
library(tidygraph)
data(interactome)
io = import('io')
sys = import('sys')

#' Return BioNet Steiner Tree subnetwork
bionet = function(g, fdr) {
    scores = setNames(pull(g, adj.p), pull(g, name))
    scores = pmax(-log10(scores) + log10(fdr), 0)
    runFastHeinz(g, scores) 
}

args = sys$cmd$parse(
    opt('g', 'gene', 'gene-level poisson', '../cis_analysis/poisson.RData'),
    opt('d', 'diff_expr', 'aneuploidy DE', '../cis_analysis/aneup_de.RData'),
    opt('o', 'outfile', 'save results .RData', 'poisson_{set}.RData'),
    opt('p', 'plotfile', 'pdf', 'bionet.pdf'))

assocs = io$load(args$gene)$result %>%
    mutate(name = toupper(external_gene_name))

net = interactome %>%
    igraph::igraph.from.graphNEL() %>%
    as_tbl_graph() %>%
    activate(nodes) %>%
    select(name = geneSymbol) %>%
    # no join, tidygraph?
    mutate(n_smp = assocs$n_smp[match(name, assocs$name)],
           n_smp = ifelse(is.na(n_smp), 0, n_smp),
           p.value = assocs$p.value[match(name, assocs$name)],
           p.value = ifelse(is.na(p.value), 1, p.value),
           adj.p = assocs$adj.p[match(name, assocs$name)],
           adj.p = ifelse(is.na(adj.p), 1, adj.p))

subnet = bionet(net, fdr=0.01)
ggraph(subnet) +
    geom_edge_link(alpha=0.2) +
    geom_node_point(aes(size=n_smp), alpha=0.2) +
    geom_node_text(aes(label = name), size=2, repel=TRUE) +
    theme_void()

dev.off()
