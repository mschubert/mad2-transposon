import_package("dplyr", attach=TRUE)
import_package("ggplot2", attach=TRUE)
import_package("BioNet", attach=TRUE)
import_package("ggraph", attach=TRUE)
import_package("tidygraph", attach=TRUE)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')

#' Return BioNet Steiner Tree subnetwork
#'
#' @param g  tidygraph-compatible network
#' @param assocs  data.frame with fields: n_smp, p.value, adj.p
#' @param thresh  p-value/fdr cutoff
#' @return  tidygraph object
bionet = function(g, assocs, thresh=0.05, var="adj.p") {
    assocs2 = assocs %>% filter(!duplicated(name))
    g = g %N>% left_join(assocs)
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

    ov2 = ov %>% select(name, p.value, statistic) %>% filter(!duplicated(name))
    net_with_stats = as_tbl_graph(net) %N>%
        left_join(ov2) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"))

    plot_net(net_with_stats, aes(size=n_smp, fill=statistic, stroke=p.value<0.05,
                 color=p.value<0.05), shape=21) +
        scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0) +
        scale_color_manual(name="signif", labels=c("n.s.", "p<0.05"),
                           values=c("white", "black"))
}

sys$run({
    args = sys$cmd$parse(
        opt('c', 'cis', 'gene-level poisson', 'poisson.rds'),
        opt('a', 'aneup', 'aneup assocs', 'ext_gene.rds'),
        opt('i', 'interactome', 'DLBCL|omnipath', 'omnipath'),
        opt('o', 'outfile', 'network data', 'bionet_omnipath.rds'),
        opt('p', 'plotfile', 'pdf', 'bionet_omnipath.pdf')
    )

    aneup = readRDS(args$aneup) %>%
        lapply(. %>% mutate(name=external_gene_name, n_smp=size, adj.p=NA))
    # Error in get.all.shortest.paths(mst, from = mst.cluster.id[j], to = mst.cluster.id[(j +  :
    #   At structural_properties.c:4863 : Weight vector must be non-negative, Invalid value
    aneup = aneup[1] # wtf...
    cis = readRDS(args$cis)$result %>% mutate(name = external_gene_name)

    if (args$interactome == "DLBCL") {
        library(DLBCL)
        data(interactome)
        net = igraph::igraph.from.graphNEL(interactome) %>%
            as_tbl_graph() %>%
            activate(nodes) %>%
            select(name = geneSymbol)
        fdr = 0.01
    } else if (args$interactome == "omnipath") {
        op = OmnipathR::import_all_interactions()
        net = op %>%
            OmnipathR::interaction_graph() %>%
            as_tbl_graph() %>%
                convert(to_undirected, .clean=TRUE) %>%
                convert(to_simple, .clean=TRUE) %>%
            mutate(name = idmap$orthologue(name, from="external_gene_name", to="mgi_symbol")) %>%
            filter(!is.na(name)) %E>%
            select(-.orig_data)
        conn = igraph::components(net)$membership == 1 # doesn't work in filter above
        net = net %N>% filter(conn) %>% activate(edges)
        fdr = 0.05
    } else
        stop("invalid interactome")

    cis_net = bionet(net, cis, fdr, "adj.p")
    ext_nets = lapply(aneup, bionet, g=net, thresh=0.1, var="p.value")
    saveRDS(list(cis_net=cis_net, ext_nets=ext_nets, aneup=aneup), file=args$outfile)

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
