library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
io = import('io')
sys = import('sys')
plt = import('plot')
idmap = import('process/idmap')

plot_pcs = function(idx, pca, x, y, hl=c()) {
    imp = summary(pca)$importance[2,c(x,y)] * 100
    pcs = names(imp)
    df = cbind(idx, pca$x) %>% mutate(ins=sample %in% hl)
    ggplot(df, aes_string(x=pcs[1], y=pcs[2])) +
        geom_point(aes(size=aneup, shape=tissue, fill=type, color=ins), stroke=1) +
        scale_shape_manual(values=c(22:30)) +
        scale_color_manual(values=c("#ffffff00", "black")) +
        geom_text_repel(aes(label=sample), color="black") +
        labs(x = sprintf("%s (%.1f%%)", pcs[1], imp[1]),
             y = sprintf("%s (%.1f%%)", pcs[2], imp[2]),
             title = "PCA plot") +
        guides(fill = guide_legend(override.aes=list(shape=21)))
}

extract_coef = function(res, coef) {
    DESeq2::lfcShrink(res, coef=coef, type="apeglm") %>%
        as.data.frame() %>%
        tibble::rownames_to_column("ensembl_gene_id") %>%
        tbl_df() %>%
        arrange(pvalue)
}

plot_volcano = function(res, coef, highlight=NULL) {
    res %>%
        mutate(label = idmap$gene(ensembl_gene_id, to="external_gene_name",
                                  dset="mmusculus_gene_ensembl"),
               circle = label %in% highlight,
               size = log10(baseMean + 1)) %>%
        plt$p_effect("padj", "log2FoldChange", thresh=0.1) %>%
        plt$volcano(p=0.1, base.size=5, label_top=30, repel=TRUE) + ggtitle(coef)
}

sys$run({
    args = sys$cmd$parse(
        opt('e', 'expr', 'gene expression RData', '../data/rnaseq/assemble.RData'),
        opt('c', 'copies', 'gene copy matrix', '../ploidy_compare/gene_copies.RData'),
        opt('a', 'aneup', 'aneuploidy score', '../ploidy_compare/analysis_set.RData'),
        opt('i', 'cis', 'cis site RData', '../cis_analysis/poisson.RData'),
        opt('o', 'outfile', 'results RData', 'aneup_de.RData'),
        opt('p', 'plotfile', 'pdf', 'aneup_de.pdf'))

    cis = io$load(args$cis)
    cis_genes = cis$result %>% filter(adj.p < 1e-3) %>% pull(external_gene_name)
    gene_copies = io$load(args$copies)
    exprset = io$load(args$expr)
    idx = exprset$idx %>%
        mutate(sample = paste0(hist_nr, tissue),
               type = sub("(^[^ ]+).*", "\\1", type)) %>%
        left_join(io$load(args$aneup)) %>%
        mutate(type = ifelse(is.na(type), "NA", type), #TODO: add annotations
               tissue = factor(tissue), type=factor(type))
    counts = exprset$counts
    narray::intersect(gene_copies, counts, along=1)
    narray::intersect(gene_copies, counts, idx$sample, along=2)

    # vst w/ copy num corr
    eset = DESeq2::DESeqDataSetFromMatrix(counts, colData=idx, ~tissue+type) %>%
        DESeq2::estimateSizeFactors(normMatrix=gene_copies)
    vs = DESeq2::getVarianceStabilizedData(DESeq2::estimateDispersions(eset))
    design(eset) = ~ tissue + type * aneup

    res = DESeq2::DESeq(eset) #, test="LRT", full=~tissue+type, reduced=~tissue)
    DESeq2::resultsNames(res)
    coefs = setdiff(DESeq2::resultsNames(res), "Intercept")
    res = sapply(coefs, extract_coef, res=res, simplify=FALSE)

    pdf(args$plotfile)
    pca = prcomp(t(vs[apply(vs, 1, var) > 0,]), center=TRUE, scale=FALSE)
    print(plot_pcs(idx, pca, 1, 2))
    print(plot_pcs(idx, pca, 3, 4))
    print(plot_pcs(idx, pca, 5, 6))
    for (name in coefs)
        print(plot_volcano(res[[name]], name, cis_genes))
    dev.off()

    save(eset, pca, cis, res, file=args$outfile)
})
