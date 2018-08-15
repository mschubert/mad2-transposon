library(cowplot)
library(patchwork)
library(dplyr)
library(magrittr)
b = import('base')
io = import('io')
seq = import('seq')
sys = import('sys')
plt = import('plot')

expr = function(segs, coords, smp, chrs=c(1:19,'X')) {
    smp = paste0("X", smp) # GRanges adds 'X' on integer start
    cur = segs %>% filter(sample == smp)
    coords = as.data.frame(coords) %>% filter(seqnames %in% chrs)
    coords$expr = coords[[smp]]

    if (nrow(cur) == 0)
        return(plot_spacer())

    ggplot(coords, aes(x=start, y=expr)) +
        geom_point(shape=1, alpha=0.3) +
        geom_hline(yintercept=1:6, color="grey", linetype="dashed") +
        geom_segment(data=cur, aes(x=start, xend=end, y=expr, yend=expr), size=3, color="green") +
        facet_grid(. ~ seqnames, scales="free_x") +
        scale_y_continuous(trans="log2", breaks=c(1:6)) +
        coord_cartesian(ylim=c(0.5,6))
}

wgs = function(segs, bins, smp, bin_y="copy.number", seg_y="mean.counts") {
    bin_ys = rlang::sym(bin_y)
    seg_ys = rlang::sym(seg_y)
    segs = segs %>% filter(Sample == smp)
    bins = bins %>% filter(Sample == smp)

    den = density(bins$counts, kernel="gaussian", bw=5)
    rpp_mode = den$x[den$y==max(den$y)] / 2

    ploidy_breaks = c(1:6) #sort(setdiff(unique(c(1,2,as.data.frame(bins)[[bin_y]])), 0))
    if (nrow(segs) == 0)
        return(plot_spacer())

    ggplot(bins) +
        geom_point(aes(x=start, y=counts), shape=1) + # x=minpoint;; start+width/2
        geom_hline(yintercept=ploidy_breaks*rpp_mode, color="grey", linetype="dashed") +
        geom_segment(data=as.data.frame(segs), size=2,
                     aes(x=start, xend=end, y=!!seg_ys, yend=!!seg_ys), color="green") +
        facet_grid(. ~ seqnames, scales="free_x") +
        scale_y_continuous(sec.axis=sec_axis(~./rpp_mode, breaks=ploidy_breaks),
                           limits=c(quantile(bins$counts, 0.01), quantile(bins$counts, 0.99)))
}

dens = function(bins, field, trans="identity", fill="blue", ...) {
    ggplot(as.data.frame(bins), aes_string(field)) +
        geom_vline(xintercept=2, linetype="dashed", alpha=0.3) +
        geom_density(fill=fill, alpha=0.5) +
        scale_x_continuous(trans=trans) +
        coord_flip(...) +
        theme(axis.title.y = element_blank())
}

plot_sample = function(smp) {
    rna_smp = tolower(sub("-(high|low)", "", smp))
    m = meta[meta$sample == rna_smp,]
    tit = sprintf("%s » %s", smp, paste(m$type, collapse=" & "))

    # DNA from 30-cell sequencing
    bins = dna$bins %>% filter(sample == smp)
    segs = dna$segments %>% filter(sample == smp)
    rpp_mode = segs$mean.counts[1] / segs$ploidy[1]
    p1 = ggplot() +
        plt$genome$pts(bins, aes(y=counts)) +
        plt$genome$segs(segs, aes(y=mean.counts), ~./rpp_mode) +
        ylab("WGS read counts") + ggtitle(tit)
    p1_dens = dens(bins %>% mutate(p = counts/rpp_mode), "p", fill="red")

    # RNA from eT ratio
    bins = rna$genes %>% filter(sample == rna_smp)
    segs = rna$segments %>% filter(sample == rna_smp)
    p2 = ggplot() +
        plt$genome$pts(bins, aes(y=expr)) +
        plt$genome$segs(segs, aes(y=expr), ~./1, breaks=1:6) +
        coord_trans(y="log2") +
        coord_cartesian(ylim=c(0.5,6)) +
        ylab("eT ratio expr")
    p2_dens = dens(bins, "expr", xlim=c(0.5,6), trans="log2")

    # Tumor weights
    pm1 = meta %>%
        mutate(sample = ifelse(meta$sample == rna_smp, rna_smp, "other")) %>%
        select(sample, Thymus=thymus_g, Spleen=spleen_g) %>%
        tidyr::gather("tissue", "grams", -sample) %>%
        mutate(category = "Weights") %>%
        ggplot(aes(x=tissue, y=grams, alpha=sample)) +
               ggbeeswarm::geom_quasirandom() +
               guides(alpha=FALSE) +
               facet_wrap(~category, scales="free_x") +
               scale_alpha_manual(values=c(1, 0.1))

    # MixCR
    pm2 = mixcr %>%
        filter(sample == rna_smp) %>%
        mutate(header = "Clonality") %>%
        ggplot(aes(x=type, y=cloneCount, fill=as.factor(cloneId))) +
               geom_col() + guides(fill=FALSE) + facet_wrap(~header)

    # Gene expression
    pm3 = expr %>%
        mutate(sample = ifelse(sample == rna_smp, rna_smp, "other"),
               category = "Gene expression") %>%
        ggplot(aes(x=gene, y=vst, alpha=sample)) +
        ggbeeswarm::geom_quasirandom() +
        facet_wrap(~category, scales="free") +
        scale_alpha_manual(values=c(1, 0.1))

    # Assemble
    plt$build_or_spacer(p1, p1_dens, p2, p2_dens, pm2, pm3)
    p = p1 + p1_dens + p2 + p2_dens +
        plot_layout(ncol=2, widths=c(10,1)) & plt$theme$no_gx()
    pm = pm1 + pm2 + pm3 +
        plot_layout(nrow=1, widths=c(2,4,length(unique(expr$gene)))) &
        theme(axis.text.x = element_text(angle=30, hjust=0.8))
    p / pm + plot_layout(heights=c(2,1))
}

if (is.null(module_name())) {
    args = sys$cmd$parse(
        opt('e', 'expr', 'gene expression RData', '../data/rnaseq/assemble.RData'),
        opt('r', 'rna', 'RNA expr reference RData', '../ploidy_from_rnaseq/eT_ploidy.RData'),
        opt('d', 'dna', '30-cell DNA WGS', '../data/wgs/30cellseq.RData'),
        opt('m', 'meta', 'sample .RData', 'analysis_set.RData'),
        opt('c', 'mixcr', 'mixcr tsv', '../data/rnaseq/mixcr_Mad2+PB.tsv'),
        opt('p', 'plotfile', 'pdf to save plot to', '/dev/null'))

    meta = io$load(args$meta)
    mixcr = io$read_table(args$mixcr, header=TRUE)
    dna = io$load(args$dna)
    rna = io$load(args$rna)

    genes = c("Mad2l1", "Msh2", "Pten", "Trp53")
    eset = io$load(args$expr)
    rownames(eset$expr) = eset$genes
    expr = narray::melt(eset$expr[genes,], dimnames=c("gene", "sample")) %>%
        mutate(gene = factor(gene, levels=genes, ordered=TRUE)) %>%
        rename(vst = value_df)

    samples = sort(unique(dna$segments$sample))
    samples = samples[sub("-.*", "", samples) %in% meta$sample]

    pdf(9, 8, file="karyograms.pdf")
    for (smp in samples) {
        message(smp)
        print(plot_sample(smp))
    }
    dev.off()
}
