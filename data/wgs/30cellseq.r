library(dplyr)
library(magrittr)
library(ggplot2)
library(patchwork)
theme_set(cowplot::theme_cowplot())
io = import('io')
seq = import('seq')
sys = import('sys')
plt = import('plot')

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
#' @return  ggplot2 object
plot_model_tracks = function(m) {
    rpp_aneufinder = as.data.frame(m$bins) %>%
        filter(copy.number != 0) %$%
        median(counts / copy.number)
    m$segments$integer_ploidy_reads = m$segments$copy.number * rpp_aneufinder
    p1 = ggplot() +
        plt$genome$pts(m$bins, aes(y=counts)) +
        plt$genome$segs(m$segments, aes(y=integer_ploidy_reads), ~./rpp_aneufinder) +
        labs(title=m$ID, subtitle="AneuFinder")

    rpp_mode = m$segments$mean.counts[1] / m$segments$ploidy[1]
    p2 = ggplot() +
        plt$genome$pts(m$bins, aes(y=counts)) +
        plt$genome$segs(m$segments, aes(y=mean.counts), ~./rpp_mode) +
        labs(subtitle="Mode diploid")

    p1 + p2 + plot_layout(ncol=1)
}

#' Plot comparison of aneuploidy scores across samples
#'
#' @param aneup  data.frame from seq$aneuploidy
#' @return  ggplot2 object
plot_aneup_compare = function(aneup) {
    aneup$sample = forcats::fct_reorder(aneup$sample, aneup$aneuploidy)

    dens = ggplot(aneup, aes(x=aneuploidy)) +
        geom_density(alpha=0.5, fill="red") +
        theme(axis.title.y=element_blank(),
              axis.text=element_blank(),
              axis.line=element_blank(),
              axis.ticks=element_blank()) +
        coord_flip()

    samples = ggplot(aneup, aes(x=sample, y=aneuploidy)) +
#        geom_segment(aes(y=sample, yend=sample, xend=aneuploidy), x=0, color="lightgrey") +
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

    mod2 = function(m, f) as.data.frame(m[[f]])
    segments = dplyr::bind_rows(lapply(models, mod2, f="segments"), .id="sample")
    bins = dplyr::bind_rows(lapply(models, mod2, f="bins"), .id="sample")
    aneup = seq$aneuploidy(segments, sample="sample", assembly="GRCm38")

    pdf(args$plotfile, 10, 6)
    print(plot_aneup_compare(aneup))
    for (m in models) {
        message(m$ID)
        print(plot_model_tracks(m))
    }
    dev.off()

    save(segments, bins, file=args$outfile)
})
