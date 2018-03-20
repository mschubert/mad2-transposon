library(ggbio)
library(dplyr)
library(magrittr)
io = import('io')

plot_cov = function(bins, segs, counts="counts", bin_y="copy.number", seg_y="mean.counts") {
    #TODO: support different name for 'counts'
    bin_ys = rlang::sym(bin_y)
    seg_ys = rlang::sym(seg_y)

    ploidy_breaks = sort(setdiff(unique(as.data.frame(bins)[[bin_y]]), 0))
    reads_per_ploidy = as.data.frame(bins) %>%
        filter(!! bin_ys != 0) %$%
        median(counts / .[[bin_y]])

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
    models = io$load("30cellseq.RData")

    pdf("non_integer.pdf", 10,4)
    for (i in seq_along(models)) {
        m = models[[i]]
        mn = sprintf("%s", m$ID) #TODO: add aneuploidy score here
        print(plot_cov(m$bins, m$segments) + ggtitle(mn))
    }
    dev.off()
})
