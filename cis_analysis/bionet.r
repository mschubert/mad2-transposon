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
    scores[is.na(scores)] = 1
    scores = pmax(-log10(scores) + log10(fdr), 0)
    as_tbl_graph(runFastHeinz(g, scores))
}

args = sys$cmd$parse(
    opt('c', 'cis', 'gene-level poisson', '../cis_analysis/poisson.RData'),
    opt('a', 'aneup', 'aneup RData', '../ploidy_compare/analysis_set.RData'),
    opt('p', 'plotfile', 'pdf', 'bionet.pdf'))

dset = io$load(args$cis)

assocs = dset$result %>%
    mutate(name = toupper(external_gene_name))

##TODO: better to plot condition-specific assoc stats on top?
#aneup = io$load(args$aneup) %>%
#    mutate(type = factor(type, levels=c("Myeloid", "T-cell", "Other"))) %>%
#    inner_join(dset$samples) %>%
#    group_by(external_gene_name) %>%
#    summarize(Myeloid = sum(type == "Myeloid" & !is.na(type)), # fraction of samples of 1 type?
#              Tcell = sum(type == "T-cell" & !is.na(type)),
#              Other = sum(type == "Other" & !is.na(type))) %>%
#    na.omit()

net = interactome %>%
    igraph::igraph.from.graphNEL() %>%
    as_tbl_graph() %>%
    activate(nodes) %>%
    select(name = geneSymbol) %>%
    # no join, tidygraph?
    mutate(n_smp = assocs$n_smp[match(name, assocs$name)],
           p.value = assocs$p.value[match(name, assocs$name)],
           adj.p = assocs$adj.p[match(name, assocs$name)])

subnet = bionet(net, fdr=0.01) #%>%
#    activate(nodes) %>%
#    mutate(Myeloid = aneup$Myeloid[match(name, toupper(aneup$external_gene_name))],
#           Tcell = aneup$Myeloid[match(name, toupper(aneup$external_gene_name))],
#           Other = aneup$Myeloid[match(name, toupper(aneup$external_gene_name))])

p = ggraph(subnet) +
    geom_edge_link(alpha=0.2) +
    geom_node_point(aes(size=n_smp), alpha=0.7) +
    geom_node_text(aes(label = name), size=2, repel=TRUE) +
    viridis::scale_color_viridis(option="magma", direction=-1) +
    theme_void()

pdf(args$plotfile)
print(p)
dev.off()
