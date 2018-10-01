library(dplyr)
library(ggplot2)
library(DESeq2)
plt = import('plot')
idmap = import('process/idmap')

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
        tibble::rownames_to_column("ensembl_gene_id") %>%
        tbl_df() %>%
        arrange(pvalue)
}

plot_volcano = function(res, highlight=NULL) {
    res %>%
        mutate(label = idmap$gene(ensembl_gene_id, to="external_gene_name",
                                  dset="mmusculus_gene_ensembl"),
               circle = label %in% highlight,
               size = log10(baseMean + 1)) %>%
        plt$p_effect("padj", "log2FoldChange", thresh=0.1) %>%
        plt$volcano(p=0.1, base.size=5, label_top=30, repel=TRUE)
}

do_wald = function(eset, fml, ex=NULL) {
    design(eset) = fml
    res = DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomWaldTest(maxit=1000)
    if (length(ex) == 0)
        ex = setdiff(DESeq2::resultsNames(res), "Intercept")
    else
        ex = grep(ex, DESeq2::resultsNames(res), value=TRUE)
    res = sapply(ex, extract_coef, res=res, simplify=FALSE)
    if (length(res) == 1)
        res = res[[1]]
    res
}

do_lrt = function(eset, fml, red) {
    design(eset) = fml
    DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomWaldTest(maxit=1000) %>%
        DESeq2::nbinomLRT(reduced=red, maxit=1000) %>%
        DESeq2::results() %>%
        as.data.frame() %>%
        tibble::rownames_to_column("ensembl_gene_id") %>%
        mutate(gene_name = idmap$gene(ensembl_gene_id, from="ensembl_gene_id",
            to="external_gene_name", dset="mmusculus_gene_ensembl")) %>%
        arrange(padj, pvalue)
}

plot_gset = function(res, sets, highlight=NULL) {
    test_one = function(set_name) {
        fdata = mutate(cur, in_set = gene_name %in% sets[[set_name]])
        mod = try(lm(stat ~ in_set, data=fdata))
        if (class(mod) == "try-error")
            return()
        broom::tidy(mod) %>%
            filter(term == "in_setTRUE") %>%
            select(-term) %>%
            mutate(size = sum(fdata$in_set, na.rm=TRUE))
    }
    cur = res %>%
        mutate(stat = log2FoldChange / lfcSE,
               gene_name = idmap$gene(ensembl_gene_id,
                    to="external_gene_name", dset="mmusculus_gene_ensembl"))
    result = sapply(names(sets), test_one, simplify=FALSE) %>%
        dplyr::bind_rows(.id="label") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
    result %>%
        plt$p_effect("adj.p", thresh=0.1) %>%
        plt$volcano(p=0.1, base.size=0.1, label_top=30, repel=TRUE, text.size=2)
}
