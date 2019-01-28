library(dplyr)
library(corrplot)
library(cowplot)
library(ggraph)
io = import('io')
sys = import('sys')
st = import('stats')

plot_cor_matrix = function(mat, title="", text_color="black") {
    p.mat = st$cor$test(mat)
    cmat = cor(mat)

    col = colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
    corrplot(cmat, method="color", col=col(200), title=title,
             order="hclust", mar=c(0,0,2,0), # title cut off otherwise
             addCoef.col = text_color, # add coefficient of correlation
             tl.col="black", tl.srt=45, #text label color and rotation
             p.mat = p.mat, sig.level = 0.05, insig = "pch", pch.col = "white"
    )
}

pcor = function(mat, fdr=1) {
    pm = GeneNet::ggm.estimate.pcor(mat, lambda=0)
    pm = GeneNet::network.test.edges(pm, direct=FALSE, plot=FALSE, verbose=TRUE)
    pm$node1 = factor(colnames(mat)[pm$node1])
    pm$node2 = factor(colnames(mat)[pm$node2])

    pm %>%
        filter(qval < fdr + .Machine$double.eps) %>%
        select(node1, node2, pcor, pval, qval) %>%
        mutate(dir=factor(sign(pcor)),
               lab=sprintf("pcor %.2f\nFDR %.2g", pcor, qval))
}

plot_pcor_net = function(pm, fdr=0.2, node_size=6, edge_size=2.5, excl=c(),
        title=sprintf("original data, FDR cutoff %.2g", fdr)) {
    g = pm %>%
        dplyr::filter(! (node1 %in% excl | node2 %in% excl)) %>%
        tidygraph::as_tbl_graph() %>%
        tidygraph::activate(edges) %>%
        tidygraph::filter(qval < fdr)

    p = ggraph(g) # no edges produce plotting error if geom_edge_link set
    if (g %>% tidygraph::as_tibble() %>% nrow() > 0)
        p = p +
            geom_edge_link(aes(label=lab, color=dir, width=abs(pcor)/10, alpha=1-qval),
                           angle_calc='along', size=edge_size)
    p = p +
        geom_node_text(aes(label=name), size=node_size) +
        theme_void() +
        ggtitle(title)

    print(p)
}

plot_bootstrapped_pcor = function(mat, fdr=0.2, n=100, show_edge_if=10, node_size=6, excl=c(),
        title=sprintf("%i bootstraps, edges if fdr<%.2f in at least %i runs", n, fdr, show_edge_if)) {
    do_bs = function(mat) {
        mat = mat[sample(seq_len(nrow(mat)), replace=TRUE),]
        pm = pcor(mat, fdr=fdr)
    }
    g = replicate(n, do_bs(mat), simplify=FALSE) %>%
        dplyr::bind_rows() %>%
        dplyr::filter(! (node1 %in% excl | node2 %in% excl)) %>%
        group_by(node1, node2) %>%
        summarize(pcor = median(pcor),
                  dir = as.factor(sign(median(pcor))),
                  n = n()) %>%
        filter(n >= show_edge_if) %>%
        tidygraph::as_tbl_graph()

    p = ggraph(g) +
        geom_edge_link(aes(label=n, color=dir, alpha=abs(pcor), width=n/10)) +#,
                       #angle_calc='along', size=2.5) +
        geom_node_text(aes(label=name), size=node_size) +
        theme_void() +
        ggtitle(title)

    print(p)
}

plot_pcor_table = function(pm, field="aneuploidy") {
    res = pm %>%
        mutate(full = paste0(node1, node2)) %>%
        filter(grepl(field, full)) %>%
        arrange(qval, pval) %>%
        transmute(node = sub(field, "", full),
                  pcor = sprintf("%.2f", pcor),
                  pval = sprintf("%.2g", pval),
                  #pval = metap::two2one(pval, invert=pcor<0),
                  fdr = sprintf("%.2g", p.adjust(pval, method="fdr")))

    gridExtra::grid.arrange(top="Full-rank partial correlations with aneuploidy (two-sided test)",
                            gridExtra::tableGrob(res))
}

sys$run({
    args = sys$cmd$parse(
        opt('e', 'expr', 'expr RData', '../expr_sets/genes_mile.RData'),
        opt('p', 'plotfile', 'pdf', 'genenet-genes_mile.pdf')
    )

    dset = io$load(args$expr)
    mat = dset$expr
    tmat = sapply(colnames(dset$groups), function(t) {
        dset$expr[,!is.na(dset$groups[,t]) & dset$groups[,t] == 1]
    }, simplify=FALSE)
    mat2 = rbind(mat, t(dset$groups[,1:3]))

    pdf(args$plotfile, 20, 15)
    plot_cor_matrix(t(mat2), text_color=NULL)
    pcor(t(mat)) %>% plot_pcor_net(node_size=4, edge_size=1)
    try(plot_bootstrapped_pcor(t(mat), node_size=4))

    pcor(t(mat2)) %>% plot_pcor_net(node_size=4, edge_size=1)
    try(plot_bootstrapped_pcor(t(mat2), node_size=4))

    for (i in seq_along(tmat)) {
        name = names(tmat)[i]
        plot_cor_matrix(t(tmat[[i]]), title=name, text_color=NULL)
        try(plot_bootstrapped_pcor(t(tmat[[i]]), fdr=0.3, node_size=4, title=name))
    }
    dev.off()
})
