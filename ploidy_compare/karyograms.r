library(ggplot2)
library(patchwork)
library(dplyr)
library(magrittr)
theme_set(cowplot::theme_cowplot())
b = import('base')
io = import('io')
seq = import('seq')
sys = import('sys')
plt = import('plot')

dens = function(bins, field, trans="identity", fill="blue", ...) {
    ggplot(as.data.frame(bins), aes_string(field)) +
        geom_vline(xintercept=2, linetype="dashed", alpha=0.3) +
        geom_density(fill=fill, alpha=0.5) +
        scale_y_continuous(trans=trans) +
        coord_flip(...) +
        theme(axis.title.y = element_blank())
}

plot_sample = function(smp) {
    rna_smp = tolower(sub("-(high|low)", "", smp))
    m = meta[meta$sample == rna_smp,]
    tit = sprintf("%s » %s", smp, paste(m$type, collapse=" & "))

    # DNA from 30-cell sequencing
    dbins = dna$bins %>% filter(sample == smp)
    dsegs = dna$segments %>% filter(sample == smp)
    rpp_mode = dsegs$mean.counts[1] / dsegs$ploidy[1]
    p1 = ggplot() +
        plt$genome$pts(dbins, aes(y=counts)) +
        plt$genome$segs(dsegs, aes(y=mean.counts), ~./rpp_mode) +
        ylab("WGS read counts") + ggtitle(tit)
    p1_dens = dens(dbins %>% mutate(p = counts/rpp_mode), "p", fill="red")

    # RNA from eT ratio
    rbins = rna$ratio %>% filter(sample == rna_smp)
    rsegs = rna$segments %>% filter(sample == rna_smp) %>% mutate(rp = ploidy/2)
    p2 = ggplot() +
        plt$genome$pts(rbins, aes(y=ratio)) +
        plt$genome$segs(rsegs, aes(y=rp), ~./0.5, breaks=1:6) +
        coord_cartesian(ylim=c(0.25,3)) +
        ylab("eT ratio expr")
    p2_dens = dens(rbins %>% mutate(p = ratio*2), "p", xlim=c(0.25,3)*2)

    # Tumor weights
    pm1 = meta %>%
        mutate(sample = ifelse(meta$sample == rna_smp, rna_smp, "other")) %>%
        select(sample, Thymus=thymus_g, Spleen=spleen_g) %>%
        tidyr::gather("tissue", "grams", -sample) %>%
        mutate(category = "Weights") %>%
        ggplot(aes(x=tissue, y=grams, alpha=sample)) +
               ggbeeswarm::geom_quasirandom(size=1.5) +
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
               sample = relevel(factor(sample), "other"),
               category = "Gene expression") %>%
        ggplot(aes(x=gene, y=counts+1, alpha=sample)) +
        geom_violin(aes(color=type), size=0.2, alpha=0.1, fill="white", position="identity", scale="width") +
        ggbeeswarm::geom_quasirandom(aes(shape=genotype), color="black", size=1.5) +
        scale_y_log10() +
        facet_wrap(~category, scales="free") +
        scale_alpha_manual(values=c(0.2, 1))

    # Assemble
    plt$build_or_spacer(p1, p1_dens, p2, p2_dens, pm2, pm3)
    p = p1 + p1_dens + p2 + p2_dens +
        plot_layout(ncol=2, widths=c(10,1)) & plt$theme$no_gx()
    widths = c(1+2, 1+length(unique(mixcr$type)), 3+length(unique(expr$gene)))
    pm = pm1 + pm2 + pm3 +
        plot_layout(nrow=1, widths=widths) &
        theme(axis.text.x = element_text(angle=40, hjust=0.8),
              axis.title.x = element_blank())
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

    meta = readRDS(args$meta)
    mixcr = io$read_table(args$mixcr, header=TRUE)
    dna = readRDS(args$dna)
    rna = readRDS(args$rna)

    genes = c("Mad2l1", "Trp53", "Ets1", "Erg", "Myc", "Sox4", "Il7", "Il7r", "Kit", "Irf4", "Stat1", "Stat3",
              "Ebf1", "Cd19", "Ighm", "Ikzf1", "Tcf15", "Gata2", "Cd55b", "Pax5", "Tmem184a", "Dntt", "Igll1",
              "Vpreb1", "Rag1", "Rag2")
    eset = io$load(args$expr)
    expr = DESeq2::DESeqDataSetFromMatrix(eset$counts, data.frame(id=colnames(eset$expr)), ~1) %>%
        DESeq2::estimateSizeFactors() %>%
        DESeq2::counts(normalized=TRUE)
    rownames(expr) = eset$genes
    expr = narray::melt(expr[genes,], dimnames=c("gene", "sample")) %>%
        mutate(gene = factor(gene, levels=genes, ordered=TRUE)) %>%
        rename(counts = value) %>%
        left_join(meta %>% select(sample, genotype, type, aneuploidy))

    samples = sort(unique(dna$segments$sample))
    samples = samples[sub("-.*", "", samples) %in% meta$sample]

    pdf(13, 10, file="karyograms.pdf")
    for (smp in samples) {
        message(smp)
        print(plot_sample(smp))
    }
    dev.off()
}
