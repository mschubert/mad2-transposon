library(ggbio)
library(dplyr)
library(magrittr)
io = import('io')
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
    # load the models we constructed
    models = io$load("../data/wgs/30cellseq.RData")
    models = models[order(names(models))]

    pdf("non_integer.pdf", 10,6)
    for (i in seq_along(models)) {
        m = models[[i]]
        message(m$ID)

        rpp_aneufinder = as.data.frame(m$bins) %>%
            filter(copy.number != 0) %$%
            median(counts / copy.number)
        p1 = plot_cov(m$bins, m$segments, reads_per_ploidy=rpp_aneufinder)

        den = density(m$bins$counts, kernel="gaussian", bw=5)
        rpp_mode = den$x[den$y==max(den$y)] / 2
        p2 = plot_cov(m$bins, m$segments, reads_per_ploidy=rpp_mode)

        p = tracks(title = sprintf("%s", m$ID), #TODO: add aneuploidy score here
                   Aneufinder = p1,
                   `Mode diploid` = p2)
        print(p)
    }
    dev.off()
})
