library(dplyr)
library(ggplot2)
library(patchwork)
seq = import('seq')
sys = import('sys')
util = import('../../ploidy_from_rnaseq/eT_ploidy')
plt = import('plot')
putil = import('../../ploidy_compare/karyograms')

# from ../../ploidy_from_rnaseq/eT_ploidy.r, but uses median instead density
extract_segment = function(smp, chr, ratio, genes) {
    `%>%` = magrittr::`%>%`
    message(smp, "", chr)
    chr_genes = genes %>%
        dplyr::filter(seqnames == chr) %>%
        dplyr::arrange(start)
    mat = ratio[chr_genes$ensembl_gene_id, smp]
    ediv = ecp::e.divisive(as.matrix(mat), min.size=50)

    res = cbind(chr_genes, cmat=mat, cluster=ediv$cluster) %>%
        dplyr::group_by(cluster) %>%
        dplyr::summarize(start = min(start),
                         end = max(end),
                         width = abs(end - start),
                         ploidy = 2 * median(cmat, na.rm=TRUE)) %>%
        dplyr::select(-cluster)
}

plot_sample = function(smp, title="") {
    rbins = cbind(genes, ratio=ratio[,smp])
    rsegs = segments %>% filter(sample == smp) %>% mutate(rp = ploidy/2)
    p2 = ggplot() +
        plt$genome$pts(rbins, aes(y=ratio)) +
        plt$genome$segs(rsegs, aes(y=rp), ~./0.5, breaks=1:6) +
        coord_cartesian(ylim=c(0.25,3)) +
        ylab("AML/normal karyotype ratio expr") +
        ggtitle(title)
    p2_dens = putil$dens(rbins %>% mutate(p = ratio*2), "p", xlim=c(0.25,3)*2)
    p2 + p2_dens + plot_layout(nrow=1, widths=c(10,1))
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'MILE study', '../../data/arrayexpress/E-GEOD-13159.rds'),
        opt('r', 'ref', 'type for diploid', 'Non-leukemia and healthy bone marrow'),
        opt('o', 'outfile', 'RData', 'aneup_scores_mile.rds'),
        opt('p', 'plotfile', 'pdf', 'aneup_scores_mile.pdf')
    )

    eset = readRDS(args$infile)
    meta = Biobase::pData(eset) %>%
        transmute(sample = rownames(.),
                  type = FactorValue..LEUKEMIA.CLASS.)
    expr = Biobase::exprs(eset)

    genes = seq$coords$gene("ensembl_gene_id", granges=TRUE) %>%
        plyranges::select(ensembl_gene_id) %>%
        as.data.frame() %>%
        filter(seqnames %in% c(1:22, 'X'))

    ref = expr[,meta$type == args$ref]
    keep = narray::map(ref, along=2, function(x) sum(x>5 & x<11) > 0.8 * length(x))
    ref = ref[keep,]
    narray::intersect(expr, ref, genes$ensembl_gene_id, along=1)
    ratio = 2^(expr - rowMeans(ref))

    # ca. 10 minutes for 2000 samples @ 100 jobs
    segments = expand.grid(sample=colnames(ratio), seqnames=c(1:22,'X'), stringsAsFactors=FALSE) %>%
        mutate(result = clustermq::Q(extract_segment, smp=sample, chr=seqnames,
                    const = list(genes=genes, ratio=ratio), n_jobs=50)) %>%
        tidyr::unnest() %>%
        group_by(sample) %>%
        mutate(ploidy = 2 + util$center_segment_density(ploidy, w=width, bw=0.5)) %>%
        ungroup()

    aneup = seq$aneuploidy(segments, sample="sample") %>%
        inner_join(meta) %>%
        arrange(-aneuploidy)

    plot_aneup = aneup %>%
        arrange(type, -aneuploidy) %>%
        group_by(type) %>%
        filter(row_number() %in% c(1:2, n()-1, n())) %>%
        ungroup()

    pdf(args$plotfile, 9, 4)
    for (i in seq_len(nrow(plot_aneup))) {
        cur = plot_aneup[i,]
        message(cur$sample)
        tit = with(cur, sprintf("%s: %s (aneup %.2f, coverage %.2f)",
                                type, sample, aneuploidy, coverage))
        print(plot_sample(cur$sample, tit) & plt$theme$no_gx())
    }
    dev.off()

    saveRDS(list(segments=segments, aneup=aneup, ratio=ratio), file=args$outfile)
})
