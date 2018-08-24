library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
io = import('io')
sys = import('sys')
plt = import('plot')
idmap = import('process/idmap')
gset = import('data/genesets')

plot_pcs = function(idx, pca, x, y, hl=c()) {
    imp = summary(pca)$importance[2,c(x,y)] * 100
    pcs = names(imp)
    df = cbind(idx, pca$x) %>% mutate(ins=sample %in% hl)
    ggplot(df, aes_string(x=pcs[1], y=pcs[2])) +
        geom_point(aes(size=aneuploidy, shape=tissue, fill=type, color=ins), stroke=1) +
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

plot_volcano = function(res, highlight=NULL) {
    res %>%
        mutate(label = idmap$gene(ensembl_gene_id, to="external_gene_name",
                                  dset="mmusculus_gene_ensembl"),
               circle = label %in% highlight,
               size = log10(baseMean + 1)) %>%
        plt$p_effect("padj", "log2FoldChange", thresh=0.1) %>%
        plt$volcano(p=0.1, base.size=5, label_top=30, repel=TRUE)
}

sys$run({
    args = sys$cmd$parse(
        opt('e', 'expr', 'gene expression RData', '../data/rnaseq/assemble.RData'),
        opt('c', 'copies', 'gene copy matrix', '../ploidy_compare/gene_copies.RData'),
        opt('m', 'meta', 'aneuploidy score', '../ploidy_compare/analysis_set.RData'),
        opt('i', 'cis', 'cis site RData', '../cis_analysis/poisson.RData'),
        opt('o', 'outfile', 'results RData', 'de.RData'),
        opt('p', 'plotfile', 'pdf', 'de.pdf'))

    cis = io$load(args$cis)
    cis_genes = cis$result %>% filter(adj.p < 1e-3) %>% pull(external_gene_name)
    gene_copies = io$load(args$copies)
    exprset = io$load(args$expr)
    idx = io$load(args$meta) %>%
        mutate(type = ifelse(is.na(type), "unknown", type))
    colnames(idx) = make.names(colnames(idx)) # fix "T-cell" in meta?
    counts = exprset$counts
    narray::intersect(gene_copies, counts, along=1)
    narray::intersect(gene_copies, counts, idx$sample, along=2)

    # vst w/ copy num corr
    eset = DESeq2::DESeqDataSetFromMatrix(counts, colData=idx, ~tissue+type) %>%
        DESeq2::estimateSizeFactors(normMatrix=gene_copies)
    vs = DESeq2::getVarianceStabilizedData(DESeq2::estimateDispersions(eset))

    # fit tissue of origin and pan-aneuploidy
    design(eset) = ~ tissue + type + aneuploidy
    robj = DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomWaldTest(maxit=1000)
    coefs = grep("tissue|aneuploidy", DESeq2::resultsNames(robj), value=TRUE)
    res = sapply(coefs, extract_coef, res=robj, simplify=FALSE)

    # fit cancer type vs mean
    cancer_type = function(term) {
        design(eset) = formula(paste("~ tissue + ", term))
        DESeq2::estimateDispersions(eset) %>%
            DESeq2::nbinomWaldTest(maxit=1000) %>%
            extract_coef(term)
    }
    ats = c("T.cell", "Myeloid", "Other") # fix "T-cell" in meta?
    res = c(res, sapply(ats, cancer_type, simplify=FALSE))

    # fit aneuploidy within cancer types
    aneup_tissue = function(type) {
        eset = eset[,eset[[type]] == 1]
        if (length(unique(eset$tissue)) == 1)
            design(eset) = formula(paste("~ aneuploidy"))
        else
            design(eset) = formula(paste("~ aneuploidy + tissue"))
        DESeq2::estimateDispersions(eset) %>%
            DESeq2::nbinomWaldTest(maxit=1000) %>%
            extract_coef("aneuploidy")
    }
    res = c(res, setNames(lapply(ats, aneup_tissue), paste0(ats, ":aneuploidy")))

    names(res) = sub("T\\.cell", "T-cell", names(res)) # fix "T-cell" in meta?

    pdf(args$plotfile)
    pca = prcomp(t(vs[apply(vs, 1, var) > 0,]), center=TRUE, scale=FALSE)
    print(plot_pcs(idx, pca, 1, 2))
    print(plot_pcs(idx, pca, 3, 4))
    print(plot_pcs(idx, pca, 5, 6))
    for (name in names(res)) {
        message(name)
        print(plot_volcano(res[[name]], cis_genes) + ggtitle(name))
    }
    dev.off()

    save(eset, pca, cis, res, file=args$outfile)
})
