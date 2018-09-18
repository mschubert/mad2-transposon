library(dplyr)
library(ggplot2)
library(patchwork)
io = import('io')
seq = import('seq')
sys = import('sys')
util = import('../../ploidy_from_rnaseq/eT_ploidy')
plt = import('plot')
putil = import('../../ploidy_compare/karyograms')

plot_sample = function(smp) {
    rbins = cbind(genes, ratio=ratio[,smp])
    rsegs = segments %>% filter(sample == smp) %>% mutate(rp = ploidy/2)
    p2 = ggplot() +
        plt$genome$pts(rbins, aes(y=ratio)) +
        plt$genome$segs(rsegs, aes(y=rp), ~./0.5, breaks=1:6) +
        coord_cartesian(ylim=c(0.25,3)) +
        ylab("AML/normal karyotype ratio expr")
    p2_dens = putil$dens(rbins %>% mutate(p = ratio*2), "p", xlim=c(0.25,3)*2)
    p2 + p2_dens + plot_layout(nrow=1, widths=c(10,1))
}

args = sys$cmd$parse(
    opt('i', 'infile', 'MILE study', '../../data/arrayexpress/E-GEOD-13159.RData'),
    opt('o', 'outfile', 'RData', 'aneup_scores_mile.RData'),
    opt('p', 'plotfile', 'pdf', 'aneup_scores_mile.pdf'))

eset = io$load(args$infile)
meta = Biobase::pData(eset) %>%
    transmute(sample = rownames(.),
              type = FactorValue..LEUKEMIA.CLASS.)
expr = Biobase::exprs(eset)

genes = seq$coords$gene("ensembl_gene_id", granges=TRUE) %>%
    plyranges::select(ensembl_gene_id) %>%
    as.data.frame() %>%
    filter(seqnames %in% c(1:22, 'X'))

ref = expr[,meta$type == "AML with normal karyotype + other abnormalities"]
keep = narray::map(ref, along=2, function(x) sum(x>5 & x<11) > 0.8 * length(x))
ref = ref[keep,]
narray::intersect(expr, ref, genes$ensembl_gene_id, along=1)
ratio = 2^(expr - rowMeans(ref))

# ca. 10 minutes for 2000 samples @ 100 jobs
segments = expand.grid(sample=colnames(ratio), seqnames=c(1:22,'X')) %>%
    mutate(result = clustermq::Q(util$extract_segment, smp=sample, chr=seqnames,
                const = list(genes=genes, ratio=ratio), n_jobs=100)) %>%
    tidyr::unnest() %>%
    group_by(sample) %>%
    mutate(ploidy = 2 + util$center_segment_density(ploidy, w=width, bw=0.25)) %>%
    ungroup()

aneup = seq$aneuploidy(segments, sample="sample") %>%
    inner_join(meta) %>%
    arrange(-aneuploidy)

plot_aneup = aneup %>%
    group_by(type) %>%
    filter(row_number() %in% c(1:2, n()-1, n())) %>%
    ungroup()

pdf(args$plotfile, 9, 4)
for (i in seq_len(nrow(plot_aneup))) {
    cur = plot_aneup[i,]
    message(cur$sample)
    tit = with(cur, sprintf("%s: %s (aneup %.2f, coverage %.2f)",
                            type, sample, aneuploidy, coverage))
    print(plot_sample(cur$sample) + ggtitle(tit) & plt$theme$no_gx())
}
dev.off()

save(segments, aneup, file=args$outfile)
