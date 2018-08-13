library(ggbio)
library(dplyr)
library(magrittr)
library(cowplot)
library(patchwork)
io = import('io')
seq = import('seq')
sys = import('sys')

#' Fit segments with assumption that modal value is diploid
#'
#' @param m  AneuFinder model object
#' @return   AneuFinder model with additional field 'ploidy' in segments
fit_segments = function(m) {
    den = density(m$bins$counts, kernel="gaussian", bw=5)
    rpp_mode = den$x[den$y==max(den$y)] / 2
    m$segments$ploidy = m$segments$mean.counts / rpp_mode
    m
}

#' Plot comparison for Aneufinder-ploidy and model-diploid ploidy
#'
#' @param m  AneuFinder model object
#' @return  ggbio object
plot_model_tracks = function(m) {
    rpp_aneufinder = as.data.frame(m$bins) %>%
        filter(copy.number != 0) %$%
        median(counts / copy.number)
    m$segments$integer_ploidy_reads = m$segments$copy.number * rpp_aneufinder
    p1 = plot_cov(m, seg_y="integer_ploidy_reads", reads_per_ploidy=rpp_aneufinder)

    rpp_mode = m$segments$mean.counts[1] / m$segments$ploidy[1]
    p2 = plot_cov(m, reads_per_ploidy=rpp_mode)

    m$segments$Sample = "constant"
    aneup = seq$aneuploidy(m$segments)$aneuploidy
    tracks(title = sprintf("%s - aneuploidy: %.2f", m$ID, aneup),
           Aneufinder = p1,
           `Mode diploid` = p2)
}

#' Plot the read coverage per bin
#'
#' @param m  AneuFinder model object
#' @param bin_y  variable name for bin counts
#' @param seg_y  variable name for for segment height
#' @param reads_per_ploidy  numeric of how to translate read counts to ploidy
#' @return  ggplot2 object
plot_cov = function(m, bin_y="copy.number", seg_y="mean.counts", reads_per_ploidy) {
    bins = m$bins
    segs = m$segments
    bin_ys = rlang::sym(bin_y)
    seg_ys = rlang::sym(seg_y)

    ploidy_breaks = sort(setdiff(unique(c(1,2,as.data.frame(bins)[[bin_y]])), 0))

    ggplot(bins) +
        geom_point(aes(x=midpoint, y=counts), shape=1) +
        geom_hline(yintercept=ploidy_breaks*reads_per_ploidy, color="grey", linetype="dashed") +
        geom_segment(data=as.data.frame(segs), size=2,
                     aes(x=start, xend=end, y=!!seg_ys, yend=!!seg_ys), color="green") +
        facet_grid(. ~ seqnames) +
        scale_y_continuous(sec.axis=sec_axis(~./reads_per_ploidy, breaks=ploidy_breaks, name="ploidy"),
                           limits=c(quantile(bins$counts, 0.01), quantile(bins$counts, 0.99))) +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              panel.spacing = unit(0, "lines"))
}

#' Plot comparison of aneuploidy scores across samples
#'
#' @param aneup  data.frame from seq$aneuploidy
#' @return  ggplot2 object
plot_aneup_compare = function(aneup) {
    aneup$Sample = forcats::fct_reorder(aneup$Sample, aneup$aneuploidy)

    dens = ggplot(aneup, aes(x=aneuploidy)) +
        geom_density(alpha=0.5, fill="red") +
        theme(axis.title.y=element_blank(),
              axis.text=element_blank(),
              axis.line=element_blank(),
              axis.ticks=element_blank()) +
        coord_flip()

    samples = ggplot(aneup, aes(x=Sample, y=aneuploidy)) +
#        geom_segment(aes(y=Sample, yend=Sample, xend=aneuploidy), x=0, color="lightgrey") +
        geom_point() +
        theme(axis.text.x = element_text(size=6, angle=90, hjust=1))

    samples + dens + plot_layout(nrow=1, widths=c(20,1))
}

sys$run({
    args = sys$cmd$parse(
        opt('o', 'outfile', 'models with copy number segments', '30cellseq.RData'),
        opt('p', 'plotfile', 'compare karyograms to aneufinder', '30cellseq.pdf'),
        arg('infiles', 'aneufinder models', sprintf('30cellseq_batch%i.RData', 1:2), arity='*'))

    # load the models we constructed
    models = do.call(c, unname(io$load(args$infiles)))
    stopifnot(sum(duplicated(names(models))) == 0)
    models = lapply(models[order(names(models))], fit_segments)

    mod2segdf = function(m) as.data.frame(m$segments)
    segments = dplyr::bind_rows(lapply(models, mod2segdf), .id="Sample")
    aneup = seq$aneuploidy(segments, sample="Sample", assembly="GRCm38")

    pdf(args$plotfile, 10, 6)
    print(plot_aneup_compare(aneup))
    for (m in models) {
        message(m$ID)
        print(plot_model_tracks(m))
    }
    dev.off()

    save(segments, file=args$outfile)
})
