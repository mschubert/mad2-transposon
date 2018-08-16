library(BioNet)
library(DLBCL)
library(ggraph)
library(tidygraph)
data(interactome)
io = import('io')
sys = import('sys')

#' Return BioNet Steiner Tree subnetwork
#'
#' @param g  tidygraph-compatible network
#' @param fdr  fdr cutoff
#' @return  tidygraph object
bionet = function(g, fdr) {
    scores = setNames(pull(g, adj.p), pull(g, name))
    scores[is.na(scores)] = 1
    scores = pmax(-log10(scores) + log10(fdr), 0)
    as_tbl_graph(runFastHeinz(g, scores))
}

#' Plot a network
#'
#' @param net  ggraph-compatible network object
#' @param node_aes  aesthetics mapping for geom_node_point
#' @return  ggplot2 object
plot_net = function(net, node_aes) {
    ggraph(net) +
        geom_edge_link(alpha=0.2) +
        geom_node_point(node_aes, alpha=0.7) +
        geom_node_text(aes(label = name), size=2, repel=TRUE) +
        viridis::scale_color_viridis(option="magma", direction=-1) +
        theme_void()
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'cis', 'gene-level poisson', '../cis_analysis/poisson.RData'),
        opt('a', 'aneup', 'aneup assocs', 'aneup_assocs.RData'),
        opt('p', 'plotfile', 'pdf', 'bionet.pdf'))

    aneup = io$load(args$aneup)
    cis = io$load(args$cis)$result %>%
        mutate(name = toupper(external_gene_name))

    net = interactome %>%
        igraph::igraph.from.graphNEL() %>%
        as_tbl_graph() %>%
        activate(nodes) %>%
        select(name = geneSymbol) %>%
        # no join, tidygraph?
        mutate(n_smp = cis$n_smp[match(name, cis$name)],
               p.value = cis$p.value[match(name, cis$name)],
               adj.p = cis$adj.p[match(name, cis$name)])

    subnet = bionet(net, fdr=0.01)

    pdf(args$plotfile)
    print(plot_net(subnet, aes(size=n_smp)))
    dev.off()
})
