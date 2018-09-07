library(dplyr)
library(tidygraph)
library(ggplot2)
library(ggraph)
io = import('io')
sys = import('sys')
st = import('stats')
idmap = import('process/idmap')
aracne = import('tools/aracne')

#TODO: move diff/sample viper to aracne module
# (also move genenet to util)
# also: split off simple cor plots in one+other+diff in separate file

#' Differential correlation of edges in a network
#'
#' @param expr  Gene expression matrix [genes x samples]
#' @param net  Coexpression network [data.frame w/ Regulator, Target]
#' @param condition  Logical vector [samples] indicating group=T/F/NA
#' @return  network with edge attributes 'diff_cor', 'diff_pval'
diff_cor = function(expr, net, condition, fdr=0.2) {
    keep = intersect(rownames(expr), net$Regulator)
    x = expr[keep ,!is.na(condition) & condition, drop=FALSE]
    y = expr[keep ,!is.na(condition) & !condition, drop=FALSE]
    de_test = function(n) wilcox.test(x[n,], y[n,]) %>% broom::tidy()
    diff_expr = sapply(rownames(x), de_test, simplify=FALSE) %>%
        dplyr::bind_rows(.id="name") %>%
        transmute(name=name, de_pval=p.value,
                  de_adjp = p.adjust(de_pval, method="fdr"))
    diff_cor = st$cor$diff_test(t(x), t(y), return_cor=TRUE)
#    res = do.call(narray::melt, diff_cor) # should work
    diff_cor = narray::melt(diff_cor$delta_cor) %>% rename(d_cor=value) %>%
        inner_join(narray::melt(diff_cor$p.value) %>% rename(p.value=value)) %>%
        rename(Regulator = Var1, Target=Var2) %>%
        arrange(p.value) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        filter(adj.p < fdr)
    g = net %>%
        select(Regulator, Target, MI) %>%
        left_join(diff_cor) %>%
        igraph::graph_from_data_frame() %>%
        tidygraph::as_tbl_graph() %>%
        left_join(diff_expr)
}

#' Condition-wise VIPER: TF activity differences
#'
#' @param expr  Gene expression matrix [genes x samples]
#' @param net  Coexpression network [igraph object]
#' @param condition  Logical vector [samples] indicating group=T/F/NA
#' @param nperm  Bootstrap permutatoins to perform (default: 1000)
#' @param ledge  Logical indicating whether to perform leading edge analysis (default: FALSE)
#' @param fdr  FDR cutoff to include filter network edges for
#' @return  data.frame with TF activity differences between groups
diff_viper = function(expr, net, condition, nperm=1000, ledge=FALSE, shadow=FALSE, fdr=0.2) {
    x = expr[,!is.na(condition) & condition, drop=FALSE]
    y = expr[,!is.na(condition) & !condition, drop=FALSE]
    sig = viper::rowTtest(x=x, y=y)
    sig = (qnorm(sig$p.value/2, lower.tail = FALSE) * sign(sig$statistic))[,1]
    sig = sig[!is.na(sig)]
    nullmodel = viper::ttestNull(x=x, y=y, per=nperm, repos=TRUE)
    regulon = aracne$to_regulon(net)
    mrs = viper::msviper(sig, regulon, nullmodel)
    if (ledge)
        mrs = viper::ledge(mrs)
    re = dplyr::arrange(dplyr::as_data_frame(summary(mrs, length(mrs[[1]]))), FDR)
    if (shadow) {
        shd = viper::shadow(mrs)$shadow
        # add shadow col to re
    }
    re
}

#' Sample-wise VIPER: TF activities per sample
#'
#' @param ...
#' @return
sample_viper = function() {
}

#' Plot a network of different activity, expression, and correlation
#'
#' @param vobj  TF activity difference data.frame from 'diff_viper'
#' @param net  Correlation network (which is subset to TFs)
#' @param diff_cor  Correlation differences between TFs data.frame from 'diff_cor'
#' @param highlight  Character vector of node identifiers to highlight
plot_subnet = function(vobj, net, fdr=0.2, highlight=c()) {
    g = net %>%
        tidygraph::as_tbl_graph() %>%
        activate(edges) %>%
        mutate(cor_dir = factor(sign(d_cor))) %>%
        activate(nodes) %>%
        inner_join(vobj %>% rename(name=Regulon)) %>%
        mutate(CIS = factor(name %in% highlight, levels=c("FALSE", "TRUE"), ordered=TRUE)) %>%
        arrange(FDR) %>%
        to_subgraph(FDR < fdr)

    ggraph(g$subgraph) +
        geom_edge_link(aes(width=MI), alpha=0.05) +
        geom_edge_link(aes(width=d_cor, color=cor_dir), alpha=0.2) +
        geom_node_point(aes(size=Size, fill=NES, stroke=de_adjp<fdr,
                            color=de_adjp<fdr, shape=CIS), alpha=0.7) +
        geom_node_text(aes(label = name), size=2, repel=TRUE) +
        scale_fill_gradient2(low="red", mid="white", high="blue", midpoint=0) +
        scale_color_manual(name="TF_de", labels=c("n.s.", paste("FDR<",fdr)),
                           values=c("#ffffff00", "#000000ff")) +
        scale_edge_color_discrete(drop=FALSE) +
        scale_shape_manual(values=c(21, 22)) +
        theme_void()
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'sample metadata', '../ploidy_compare/analysis_set.RData'), # only mad2pb
        opt('e', 'expr', 'expr RData', '../data/rnaseq/assemble.RData'),
        opt('n', 'network', 'aracne', '../data/networks/E-GEOD-13159.RData'),
        opt('c', 'cis', 'common insertion .RData', '../cis_analysis/poisson.RData'),
        opt('f', 'fdr', 'CIS fdr for highlight', '0.001'),
        opt('p', 'plotfile', 'pdf', 'viper_mad2pb.pdf')
    )

    dset = io$load(args$expr)
    net = io$load(args$network)
    tf_net = filter(net, Target %in% Regulator)

    if (grepl("rnaseq/assemble.RData", args$expr)) {
        types = dset$idx$type
        expr = dset$expr
        rownames(expr) = idmap$orthologue(rownames(expr),
            from="ensembl_gene_id", to="hgnc_symbol", dset="mmusculus_gene_ensembl")
        expr = expr[!is.na(rownames(expr)),]
        highlight = io$load(args$cis)$result %>%
            filter(adj.p < as.numeric(args$fdr)) %>%
            pull(external_gene_name) %>%
            idmap$orthologue(from="external_gene_name", to="hgnc_symbol",
                             dset="mmusculus_gene_ensembl")
    } else {
        types = Biobase::pData(dset)$FactorValue..LEUKEMIA.CLASS.
        expr = Biobase::exprs(dset)
        highlight = c()
    }

    tmat = narray::mask(types, along=2) + 0
    if (grepl("rnaseq/assemble.RData", args$expr)) {
        meta = io$load(args$meta)
        meta = meta[match(colnames(expr), meta$sample),]
        mask = cbind(all=1, tmat)
        colnames(mask) = paste("aneuploidy", colnames(mask), sep=":")
        mask[!mask] = NA
        add = mask * meta$aneuploidy
        mask = apply(add, 2, function(x) x > median(x, na.rm=TRUE)) + 0
        tmat = cbind(tmat, mask)
    }

    pdf(args$plotfile)
    for (cond in colnames(tmat)) {
        dviper = diff_viper(expr, net, tmat[,cond])
        dcor = diff_cor(expr, tf_net, tmat[,cond])
        p = try(plot_subnet(dviper, dcor, highlight=highlight))
        if (class(p) == "try-error" || class(try(ggplot_build(p))) == "try-error")
            p = ggplot(data.frame()) + geom_point() + xlim(0, 10) + ylim(0, 100)
        print(p + ggtitle(cond))
    }
    dev.off()
})
