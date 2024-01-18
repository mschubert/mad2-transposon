library(dplyr)
library(SummarizedExperiment)
sys = import('sys')
plt = import('plot')
gset = import('genesets')
idmap = import('process/idmap')
deseq = import('process/deseq')

load_counts = function(count_file, samples) {
    seq_lib = tools::file_path_sans_ext(basename(count_file))
    message(seq_lib)
    counts = readRDS(count_file)[,samples$sample,drop=FALSE]
    colnames(counts) = sprintf("%s.%s", seq_lib, samples$sample)
    counts
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

sys$run({
    args = sys$cmd$parse(
        opt('s', 'samples', 'tsv', 'samples.tsv'),
        opt('o', 'outfile', 'rds', 'eset.rds'),
        opt('p', 'plotfile', 'pdf', 'eset.pdf'),
        arg('infiles', 'rds', arity='*',
            list.files('seq_counts', '\\.rds$', recursive=TRUE, full.names=TRUE))
    )

    meta = readr::read_tsv(args$samples, comment="#") |>
        mutate(short = sprintf("%s %s %s %sh-%i", cline, genotype, conc, hours, rep),
               rep = factor(rep))
    counts = lapply(args$infiles, load_counts, samples=meta) |>
        narray::stack(along=2, fill=0)

    eset = DESeq2::DESeqDataSetFromMatrix(counts, meta, ~1)
    vs = DESeq2::varianceStabilizingTransformation(eset)

    pdf(args$plotfile, 12, 10)
    print(plot_pca(vs, Inf))
    print(plot_pca(vs, 500))
    print(plot_hallmark_gsva(vs))
    dev.off()

    saveRDS(eset, file=args$outfile)
})
