import_package("dplyr", attach=TRUE)
import_package("ggplot2", attach=TRUE)
import_package("DESeq2", attach=TRUE)
gset = import('genesets')
plt = import('plot')

plot_pcs = function(idx, pca, x, y, hl=c()) {
    imp = summary(pca)$importance[2,c(x,y)] * 100
    pcs = names(imp)
    df = cbind(idx, pca$x) %>% mutate(ins=sample %in% hl)
    ggplot(df, aes_string(x=pcs[1], y=pcs[2])) +
        geom_point(aes(size=aneuploidy, shape=tissue, fill=type, color=ins), stroke=1) +
        scale_shape_manual(values=c(22:30)) +
        scale_color_manual(values=c("#ffffff00", "black")) +
        ggrepel::geom_text_repel(aes(label=sample), color="black") +
        labs(x = sprintf("%s (%.1f%%)", pcs[1], imp[1]),
             y = sprintf("%s (%.1f%%)", pcs[2], imp[2]),
             title = "PCA plot") +
        guides(fill = guide_legend(override.aes=list(shape=21)))
}

extract_coef = function(res, coef, type="apeglm") {
    DESeq2::lfcShrink(res, coef=coef, type=type) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene_name") %>%
        as_tibble() %>%
        mutate(stat = log2FoldChange / lfcSE) %>%
        arrange(padj, pvalue)
}

plot_volcano = function(res, highlight=NULL) {
    if ("baseMean" %in% colnames(res))
        df = res %>%
            mutate(label = gene_name,
                   circle = label %in% highlight,
                   size = log10(baseMean + 1)) %>%
            plt$p_effect("padj", "log2FoldChange", thresh=0.1)
    else
        df = res %>%
            mutate(label = gene_name,
                   circle = label %in% highlight,
                   size = mean_expr) %>%
            plt$p_effect("adj.p", "estimate", thresh=0.1)

    plt$volcano(df, p=0.1, base.size=5, label_top=30, repel=TRUE)
}

do_wald = function(eset, fml, ex=NULL, drop=TRUE) {
    design(eset) = fml
    res = DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomWaldTest(maxit=1000)
    if (length(ex) == 0)
        ex = setdiff(DESeq2::resultsNames(res), "Intercept")
    else
        ex = grep(ex, DESeq2::resultsNames(res), value=TRUE)
    res = sapply(ex, extract_coef, res=res, simplify=FALSE)
    if (length(res) == 1 && drop)
        res = res[[1]]
    res
}

do_lrt = function(eset, fml, red) {
    design(eset) = fml
    DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomLRT(reduced=red, maxit=1000) %>%
        DESeq2::results() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene_name") %>%
        arrange(padj, pvalue)
}

plot_gset = function(res, sets, highlight=NULL, fdr=0.1, base.size=0.1,
                     label_top=30, repel=TRUE, stat="stat") {
    #todo: it seems better to do gset$test_lm with stat=log2FoldChange if the
    # FCs were shrunk before; otherwise, for a Wald test stat=stat looks better
    # and for an LRT sign(lfc)*stat
    p = res %>%
        gset$test_lm(sets, stat=stat) %>%
        plt$p_effect("adj.p", thresh=fdr) %>%
        plt$volcano(p=fdr, base.size=base.size, label_top=label_top,
                    repel=repel, text.size=2)
    plt$try(p)
}
