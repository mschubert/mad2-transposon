library(dplyr)
library(SummarizedExperiment)
sys = import('sys')
plt = import('plot')
gset = import('genesets')
idmap = import('process/idmap')
deseq = import('process/deseq')

read_counts = function(rfiles, sel="sense") {
    samples = sub(".ReadsPerGene.out.tab", "", basename(rfiles), fixed=TRUE)
    reads = lapply(rfiles, readr::read_tsv, col_names=FALSE) %>%
        setNames(samples) %>%
        bind_rows(.id="sample") %>%
        dplyr::rename(feature=X1, unstranded=X2, sense=X3, antisense=X4)

    counts = reads %>% filter(!grepl("^N_", feature)) %>%
        select(sample, feature, {{ sel }}) %>%
        tidyr::pivot_wider(names_from=sample, values_from={{ sel }})
    cmat = data.matrix(counts[-1])
    rownames(cmat) = counts$feature
    cmat
}

plot_pca = function(eset, ntop=500, title="PCA") {
    pr = deseq$prcomp(eset, ntop=ntop)
    lines = as_tibble(as.data.frame(colData(eset))) |>
        mutate(PC1 = pr$x[,"PC1"], PC2=pr$x[,"PC2"])
    lines2 = lines |>
        group_by(sample) |>
        filter(n_distinct(rep) > 1)

    plt$pca(eset, aes(x=PC1, y=PC2, fill=paste(cline, genotype), color=rep), pr=pr) +
        geom_line(data=lines, aes(group=treatment), color="black", alpha=0.1) +
        geom_point(aes(shape=factor(conc)), size=3, color="#12121256") +
        scale_shape_manual(values=c("0"=21, "250"=22, "500"=23)) +
        scale_fill_brewer(palette="Paired") +
        scale_color_brewer(palette="Dark2") +
        guides(fill = guide_legend(override.aes=list(shape=21))) +
        ggrepel::geom_text_repel(aes(label=short)) +
        labs(title = title,
             subtitle = sprintf("%i most variable genes", length(pr$center)))
}

plot_hallmark_gsva = function(vs) {
    sets = gset$get_human("MSigDB_Hallmark_2020")
    rownames(vs) = idmap$gene(rownames(vs), to="hgnc_symbol")
    colnames(vs) = vs$short
    vs = vs[!is.na(rownames(vs)),]
    scores = GSVA::gsva(assay(vs), sets)
    ComplexHeatmap::Heatmap(scores, row_dend_reorder=TRUE)
}

plot_qPCR_genes = function(vs) {
    rownames(vs) = idmap$gene(rownames(vs), to="hgnc_symbol")
    colnames(vs) = vs$short
    dset = vs[c("CXCL10", "ISG15", "IFNB1", "IL6", "CXCL8", "CCL5"),] |>
        DESeq2::counts() |> t() |> cbind(colData(vs)) |>
        as.data.frame() |> as_tibble() |>
        tidyr::pivot_longer(-(sample:short), names_to="gene", values_to="counts")
    ggplot(dset, aes(x=short, y=counts)) +
        geom_col(aes(fill=factor(conc))) +
        scale_y_continuous(trans="log1p", breaks=c(1,2,5,20,100,500)) +
        facet_wrap(~ gene) +
        theme(axis.text.x = element_text(angle=90, hjust=1))
}

sys$run({
    args = sys$cmd$parse(
        opt('s', 'samples', 'tsv', 'FF230302.tsv'),
        opt('o', 'outfile', 'rds', 'FF230302.rds'),
        opt('p', 'plotfile', 'pdf', 'FF230302.pdf'),
        arg('infiles', '.ReadsPerGene.out.tab', arity='*',
            list.files("seq_aligned/FF230302", "\\.ReadsPerGene\\.out\\.tab", full.names=TRUE))
    )

    meta = readr::read_tsv(args$samples, comment="#") |>
        mutate(short = sprintf("%s %s %s %sh-%i", cline, genotype, conc, hours, rep),
               rep = factor(rep))
    counts = read_counts(args$infiles)[,meta$sample]
    eset = DESeq2::DESeqDataSetFromMatrix(counts, meta, ~1)
    vs = DESeq2::varianceStabilizingTransformation(eset)

    pdf(args$plotfile, 12, 10)
    print(plot_pca(vs, Inf))
    print(plot_pca(vs, 500))
    print(plot_hallmark_gsva(vs))
    print(plot_hallmark_gsva(eset[,vs$cline == "BT549" & vs$genotype == "WT"]))
    print(plot_qPCR_genes(eset[,eset$cline == "BT549" & eset$genotype == "WT"]))
    dev.off()

    saveRDS(eset, file=args$outfile)
})
