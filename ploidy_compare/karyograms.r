library(cowplot)
library(patchwork)
library(dplyr)
library(magrittr)
b = import('base')
io = import('io')
seq = import('seq')
sys = import('sys')
plt = import('plot')

dens = function(bins, field, trans="identity", fill="blue", ...) {
    ggplot(as.data.frame(bins), aes_string(field)) +
        geom_vline(xintercept=2, linetype="dashed", alpha=0.3) +
        geom_density(fill=fill, alpha=0.5, bw=0.5) +
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
               sample = relevel(factor(sample), "other"),
               category = "Gene expression") %>%
        ggplot(aes(x=gene, y=vst, alpha=sample)) +
        ggbeeswarm::geom_quasirandom() +
        facet_wrap(~category, scales="free") +
        scale_alpha_manual(values=c(0.1, 1))

    # Assemble
    plt$build_or_spacer(p1, p1_dens, p2, p2_dens, pm2, pm3)
    p = p1 + p1_dens + p2 + p2_dens +
        plot_layout(ncol=2, widths=c(10,1)) & plt$theme$no_gx()
    pm = pm1 + pm2 + pm3 +
        plot_layout(nrow=1, widths=c(2,4,length(unique(expr$gene)))) &
        theme(axis.text.x = element_text(angle=30, hjust=0.8),
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
