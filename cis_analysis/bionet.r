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
#' @param assocs  data.frame with fields: n_smp, p.value, adj.p
#' @param thresh  p-value/fdr cutoff
#' @return  tidygraph object
bionet = function(g, assocs, thresh=0.05, var="adj.p") {
    g = g %>% activate(nodes) %>% left_join(assocs)
    scores = setNames(pull(g, !! rlang::sym(var)), pull(g, name))
    scores[is.na(scores)] = 1
    scores = pmax(-log10(scores) + log10(thresh), 0)
    as_tbl_graph(runFastHeinz(g, scores))
}

#' Plot a network
#'
#' @param net  ggraph-compatible network object
#' @param node_aes  aesthetics mapping for geom_node_point
#' @return  ggplot2 object
plot_net = function(net, node_aes, ...) {
    set.seed(121979) # same layout if same nodes
    p = ggraph(net) +
        geom_node_point(node_aes, alpha=0.7, ...) +
        geom_node_text(aes(label = name), size=2, repel=TRUE) +
        viridis::scale_fill_viridis(option="magma", direction=-1) +
        theme_void()
    if (igraph::gsize(net) > 0)
        p = p + geom_edge_link(alpha=0.2)
    p
}

#' Plot a network overlayed with external associations
#'
#' @param net  ggraph-compatible network object
#' @param ov  data.frame with fields: external_gene_name, statistic
#' @return  ggplot2 object
plot_net_overlay = function(net, ov) {
    if (igraph::vcount(net) == 0)
        return(patchwork::plot_spacer())
    net %>%
        as_tbl_graph() %>%
        activate(nodes) %>%
        mutate(p.value = ov$p.value[match(name, toupper(ov$external_gene_name))],
               statistic = ov$statistic[match(name, toupper(ov$external_gene_name))],
               adj.p = p.adjust(p.value, method="fdr")) %>%
        plot_net(aes(size=n_smp, fill=statistic, stroke=p.value<0.05,
                     color=p.value<0.05), shape=21) +
            scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0) +
            scale_color_manual(name="signif", labels=c("n.s.", "p<0.05"),
                               values=c("white", "black"))
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'cis', 'gene-level poisson', '../cis_analysis/poisson.RData'),
        opt('a', 'aneup', 'aneup assocs', 'ext_gene.RData'),
        opt('p', 'plotfile', 'pdf', 'bionet.pdf'))

    aneup = io$load(args$aneup)
    cis = io$load(args$cis)$result %>%
        mutate(name = toupper(external_gene_name))

    net = interactome %>%
        igraph::igraph.from.graphNEL() %>%
        as_tbl_graph() %>%
        activate(nodes) %>%
        select(name = geneSymbol)

    cis_net = bionet(net, cis, 0.01, "adj.p")
    ext_nets = lapply(aneup, function(a) {
        bionet(net, mutate(a, name=toupper(external_gene_name), adj.p = NA, n_smp=size),
               thresh=0.1, var="p.value")
    })

    pdf(args$plotfile)
    print(plot_net(cis_net, aes(size=n_smp)))
    for(ov in names(aneup)) {
        message(ov)
        print(plot_net_overlay(cis_net, aneup[[ov]]) + ggtitle(ov))
        print(plot_net_overlay(ext_nets[[ov]], aneup[[ov]]) +
              ggtitle(sprintf("%s (subnetwork)", ov)))
    }
    dev.off()
})
