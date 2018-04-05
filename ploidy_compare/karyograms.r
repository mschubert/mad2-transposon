library(cowplot)
library(patchwork)
library(dplyr)
library(magrittr)
b = import('base')
io = import('io')
seq = import('seq')

#args = sys$cmd$parse(
#    opt('i', 'infile', 'RData aneuploidy scores', 'copy_segments.RData'),
#    opt('p', 'plotfile', 'pdf to save plot to', '/dev/null'))

dna = io$load("../ploidy_from_wgs/copy_segments.RData")
dna_aneup = dna$segments %>%
    seq$aneuploidy() %>%
    dplyr::rename(sample = Sample) %>%
    mutate(sample = sub("(-high)|(-low)", "", sample),
           sample = paste0(sub("[^0-9]+", "", sample), sub("[^ST]+", "", sample)))
rna = io$load("../ploidy_from_rnaseq/panel_ploidy.RData")
rna_aneup = rna$segments %>%
    seq$aneuploidy(sample="sample", ploidy="expr") %>%
    mutate(sample = toupper(gsub("[^0-9stST]+", "", sample)))
eT2 = io$load("../ploidy_from_rnaseq/eT_ploidy.RData")
eT2_aneup = eT2$segments %>%
    seq$aneuploidy(sample="sample", ploidy="expr") %>%
    mutate(sample = toupper(gsub("[^0-9stST]+", "", sample)))
eT_aneup = io$load("../tis_rna-ploidy_cis-fet/aneuploidy_mad2.RData") %>%
    transmute(sample=sample, aneuploidy=aneup)

tissues = setNames(c("spleen", "thymus"), c("S","T"))

aneups = list(WGS = dna_aneup, `RNA-seq (eDivisive)` = rna_aneup,
              `RNA-seq (eT)` = eT_aneup, `RNA-seq (eT, new)` = eT2_aneup) %>%
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

    expr = function(segs, coords, smp) {
        smp = sub("-(high|low)", "", smp)
        smp = paste0(sub("[^ST]+", "", smp), sub("[^0-9]+", "", smp))
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
        reads_per_ploidy = as.data.frame(bins) %>%
            filter(copy.number != 0) %$%
            median(counts / copy.number)
        ploidy_breaks = c(1:6) #sort(setdiff(unique(c(1,2,as.data.frame(bins)[[bin_y]])), 0))
        segs = segs %>% filter(Sample == smp)
        if (nrow(segs) == 0)
            return(plot_spacer())
        ggplot(bins %>% filter(Sample == smp)) +
            geom_point(aes(x=start, y=counts), shape=1) + # x=minpoint;; start+width/2
            geom_hline(yintercept=ploidy_breaks*reads_per_ploidy, color="grey", linetype="dashed") +
            geom_segment(data=as.data.frame(segs), size=2,
                         aes(x=start, xend=end, y=!!seg_ys, yend=!!seg_ys), color="green") +
            facet_grid(. ~ seqnames, scales="free_x") +
            scale_y_continuous(sec.axis=sec_axis(~./reads_per_ploidy, breaks=ploidy_breaks, name="ploidy"),
                               limits=c(quantile(bins$counts, 0.01), quantile(bins$counts, 0.99)))
    }

    p1 = expr(eT2$segments, eT2$genes, smp) + ylab("eT ratio expr") + ggtitle(smp)
    p2 = expr(rna$segments, rna$genes, smp) + ylab("RNA panel expr")
    p3 = wgs(dna$segments, dna$bins, smp) + ylab("WGS read counts")
    p = p1 + p2 + p3 + plot_layout(ncol=1) & mytheme
}

pdf(9, 9, file="karyograms.pdf")
for (smp in unique(dna$segments$Sample)) {
    message(smp)
    print(plot_sample(smp))
}
dev.off()
