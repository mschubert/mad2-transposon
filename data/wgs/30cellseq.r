library(ggbio)
library(dplyr)
library(magrittr)
io = import('io')
seq = import('seq')
sys = import('sys')

plot_cov = function(bins, segs, counts="counts", bin_y="copy.number", seg_y="mean.counts",
                    reads_per_ploidy) {
    #TODO: support different name for 'counts'
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

sys$run({
    args = sys$cmd$parse(
        opt('o', 'outfile', 'models with copy number segments', '30cellseq.RData'),
        opt('p', 'plotfile', 'compare karyograms to aneufinder', '30cellseq.pdf'),
        arg('infiles', 'aneufinder models', sprintf('30cellseq_batch%i.RData', 1:2), arity='*'))

    # load the models we constructed
    models = do.call(c, unname(io$load(args$infiles)))
    stopifnot(sum(duplicated(names(models))) == 0)
    models = models[order(names(models))]
    bins = lapply(models, function(m) as.data.frame(m$bins)) %>%
        dplyr::bind_rows(.id="Sample") %>%
        tbl_df()
    segments = list()

    pdf(args$plotfile, 10,6)
    for (i in seq_along(models)) {
        m = models[[i]]
        message(m$ID)

        rpp_aneufinder = as.data.frame(m$bins) %>%
            filter(copy.number != 0) %$%
            median(counts / copy.number)
        m$segments$integer_ploidy_reads = m$segments$copy.number * rpp_aneufinder
        p1 = plot_cov(m$bins, m$segments, seg_y="integer_ploidy_reads", reads_per_ploidy=rpp_aneufinder)

        den = density(m$bins$counts, kernel="gaussian", bw=5)
        rpp_mode = den$x[den$y==max(den$y)] / 2
        p2 = plot_cov(m$bins, m$segments, reads_per_ploidy=rpp_mode)

        segments[[m$ID]] = as.data.frame(m$segments) %>%
            mutate(Sample = m$ID,
                   ploidy = mean.counts / rpp_mode)
        aneup = seq$aneuploidy(segments[[m$ID]])

        p = tracks(title = sprintf("%s - aneuploidy: %.2f", m$ID, aneup$aneuploidy),
                   Aneufinder = p1,
                   `Mode diploid` = p2)
        print(p)
    }
    dev.off()

    segments = do.call(bind_rows, segments)
    save(bins, segments, file=args$outfile)
})
