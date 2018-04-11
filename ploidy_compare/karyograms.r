library(cowplot)
library(patchwork)
library(dplyr)
library(magrittr)
b = import('base')
io = import('io')
seq = import('seq')
sys = import('sys')

args = sys$cmd$parse(
    opt('e', 'euploid', 'RNA expr reference RData', '../ploidy_from_rnaseq/eT_ploidy.RData'),
    opt('d', 'dna', '30-cell DNA WGS', '../ploidy_from_wgs/copy_segments.RData'),
    opt('m', 'meta', 'sample .RData', '../data/meta/meta.RData'),
    opt('p', 'plotfile', 'pdf to save plot to', '/dev/null'))

meta = io$load(args$meta)
dna = io$load(args$dna)
dna_aneup = dna$segments %>%
    seq$aneuploidy() %>%
    dplyr::rename(sample = Sample) %>%
    mutate(sample = sub("(-high)|(-low)", "", sample),
           sample = paste0(sub("[^0-9]+", "", sample), sub("[^ST]+", "", sample)))
eT = io$load(args$euploid)
eT_aneup = eT$segments %>%
    seq$aneuploidy(sample="sample", ploidy="expr") %>%
    mutate(sample = toupper(gsub("[^0-9stST]+", "", sample)))

tissues = setNames(c("spleen", "thymus"), c("S","T"))

aneups = list(WGS = dna_aneup, `RNA-seq (eT, eDivisive)` = eT_aneup) %>%
    dplyr::bind_rows(.id="type") %>%
    group_by(type) %>%
    mutate(aneuploidy = scale(aneuploidy, center=FALSE)) %>%
    ungroup() %>%
    mutate(sample = forcats::fct_reorder(sample, aneuploidy),
           tissue = tissues[sub("Healthy|[0-9]+", "", sample)])

plot_sample = function(smp, chrs=c(1:19,'X')) {
    mytheme = theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0, "lines"))
    notheme = theme(
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

    expr = function(segs, coords, smp) {
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
            notheme
    }

    hist = sub("[^0-9]+", "", smp)
    rna_smp = sub("-(high|low)", "", smp)
    rna_smp = paste0(sub("[^ST]+", "", rna_smp), hist)
    m = meta$meta[meta$meta$`Hist nr.` == hist,]
    if (nrow(m) == 1)
        tit = sprintf("%s » %s", smp, m$Diagnosis)
    else
        tit = smp

    p1 = wgs(dna$segments, dna$bins, smp) + ylab("WGS read counts") + ggtitle(tit)
    dna_smp = dna$bins %>% filter(Sample == smp)
    p1_dens = dens(dna_smp, "counts", fill="red",
            xlim=c(quantile(dna_smp$counts, 0.01), quantile(dna_smp$counts, 0.99))) +
        theme(axis.title.y = element_blank())
    p2 = expr(eT$segments, eT$genes, rna_smp) + ylab("eT ratio expr")
    p2_dens = dens(eT$genes, rna_smp, xlim=c(0.5,6), trans="log2")
    if (class(try(ggplot_build(p2_dens))) == "try-error")
        p2_dens = plot_spacer()
    p = p1 + p1_dens + p2 + p2_dens + plot_layout(ncol=2, widths=c(10,1)) & mytheme

    meta$meta$sample = ifelse(meta$meta$`Hist nr.` == hist, rna_smp, "other")
    weights = meta$meta %>%
        transmute(sample = sample,
                  type = "Tumor weight",
                  thymus = as.numeric(`Thymus (mg)`)/1000,
                  spleen = as.numeric(`Spleen (mg)`)/1000) %>%
        tidyr::gather("tissue", "grams", -sample, -type)
    icr = meta$icr %>%
        filter(id == rna_smp) %>%
        mutate(header = "Clonality")
    meta$expression$sample = ifelse(meta$expression$id == rna_smp, rna_smp, "other")
    meta$expression$type = "Gene expression"
    pm1 = ggplot(weights, aes(x=tissue, y=grams, alpha=sample)) +
        ggbeeswarm::geom_quasirandom() +
        guides(alpha=FALSE) +
        facet_wrap(~type) +
        scale_alpha_manual(values=c(0.1, 1))
#    if (sum(weights$sample != "other") == 0)
#        pm1 = plot_spacer()
    pm2 = ggplot(icr, aes(x=type, y=counts, fill=gene)) +
        geom_col() + guides(fill=FALSE) + facet_wrap(~header)
    if (class(try(ggplot_build(pm2))) == "try-error")
        pm2 = plot_spacer()
    pm3 = ggplot(meta$expression, aes(x=gene, y=vst, alpha=sample)) +
        ggbeeswarm::geom_quasirandom() +
        facet_wrap(~type, scales="free") +
        scale_alpha_manual(values=c(0.1, 1))
    if (class(try(ggplot_build(pm3))) == "try-error" || sum(meta$expression$sample != "other") == 0)
        pm3 = plot_spacer()
    pm = pm1 + pm2 + pm3 + plot_layout(nrow=1, widths=c(2,1,length(unique(meta$expression$gene))))

    p / pm + plot_layout(heights=c(2,1))
}

pdf(9, 8, file="karyograms.pdf")
for (smp in unique(dna$segments$Sample)) {
    message(smp)
    print(plot_sample(smp))
}
dev.off()
