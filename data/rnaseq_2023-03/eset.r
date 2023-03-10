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
    idx = samples[[seq_lib]]
    counts = readRDS(count_file)[,names(idx),drop=FALSE]
    colnames(counts) = sprintf("%s.%s", seq_lib, names(idx))
    counts
}

plot_pca = function(eset, ntop=500, title="PCA") {
    pr = deseq$prcomp(eset, ntop=ntop)
    lines = as_tibble(as.data.frame(colData(eset))) %>%
        mutate(PC1 = pr$x[,"PC1"], PC2=pr$x[,"PC2"])
    lines2 = lines %>%
        group_by(sample) %>%
        filter(n_distinct(batch) > 1)

    plt$pca(eset, aes(x=PC1, y=PC2, fill=short, color=batch), pr=pr) +
        geom_line(data=lines, aes(group=short), color="black", alpha=0.1) +
        geom_point(shape=21, size=3, color="#12121256") +
        scale_fill_brewer(palette="Paired") +
        scale_color_brewer(palette="Dark2") +
        ggrepel::geom_text_repel(aes(label=sample)) +
        labs(title = title,
             subtitle = sprintf("%i most variable genes", length(pr$center)))
}

plot_hallmark_gsva = function(vs) {
    sets = gset$get_human("MSigDB_Hallmark_2020")
    rownames(vs) = idmap$gene(rownames(vs), to="hgnc_symbol")
    colnames(vs) = eset$sample
    vs = vs[!is.na(rownames(vs)),]
    scores = GSVA::gsva(assay(vs), sets)
    ComplexHeatmap::Heatmap(scores, row_dend_reorder=TRUE)
}

sys$run({
    args = sys$cmd$parse(
        opt('s', 'samples', 'yaml', 'fastq_id.yaml'),
        opt('o', 'outfile', 'rds', 'eset.rds'),
        opt('p', 'plotfile', 'pdf', 'eset.pdf'),
        arg('infiles', 'rds', arity='*',
            list.files('seq_counts', '\\.rds$', recursive=TRUE, full.names=TRUE))
    )

    samples = yaml::read_yaml("fastq_id.yaml") %>% lapply(unlist)
    counts = lapply(args$infiles, load_counts, samples=samples) %>%
        narray::stack(along=2, fill=0)

    meta = tibble(sample = unlist(samples, use.names=FALSE),
                  batch = sub("^(.+)\\.[^.]+", "\\1", colnames(counts))) %>%
        mutate(genotype = sub("([^_]+)_([^-]+)-(.+)", "\\1", sample),
               treatment = sub("([^_]+)_([^-]+)-(.+)", "\\2", sample),
               rep = sub("([^_]+)_([^-]+)-(.+)", "\\3", sample),
               short = sub("([^_]+)_([^-]+)-(.+)", "\\1_\\2", sample))
    eset = DESeq2::DESeqDataSetFromMatrix(counts, meta, ~1)
    vs = DESeq2::varianceStabilizingTransformation(eset)

    pdf(args$plotfile, 12, 10)
    print(plot_pca(vs, Inf))
    print(plot_pca(vs, 500))
    print(plot_hallmark_gsva(vs))
    dev.off()

    saveRDS(eset, file=args$outfile)
})
