library(dplyr)
library(patchwork)
library(flowCore)
plt = import('plot')

ggfacs = function(ff, aes) {
    meta = ff@parameters@data
    df = as_tibble(ff@exprs)
    colnames(df) = ifelse(is.na(meta$desc), meta$name, meta$desc)
    comp = ff@description$SPILL
    name = ff@description$`TUBE NAME`

    ggplot(df, aes) +
        geom_bin2d(bins = 70) +
        scale_fill_continuous(type="viridis", trans="log") +
        theme_bw() +
        xlim(0, 2.5e5) +
        ylim(0, 2.5e5) +
        ggtitle(sprintf("%s vs. %s", all.vars(aes[[1]]), all.vars(aes[[2]])))
}

plot_one = function(fname) {
    fc = read.FCS(fname)
    bn = tools::file_path_sans_ext(basename(fname))
    bn = sub("Specimen_002_", "", bn, fixed=TRUE)
    plots = list(
        ggfacs(fc, aes(x=`FSC-H`, y=`SSC-H`)),
        ggfacs(fc, aes(x=`CD45`, y=`FSC-H`)) + scale_x_log10(),
        ggfacs(fc, aes(x=`MAC1/GR1`, y=`SSC-A`)) + scale_x_log10(),
        ggfacs(fc, aes(x=`FSC-A`, y=`B220`)) + scale_y_log10(),
        ggfacs(fc, aes(x=`CD3`, y=`FSC-H`)) + scale_x_log10(),
        ggfacs(fc, aes(x=`SCA-1`, y=`CKIT`)) + scale_x_log10() + scale_y_log10()
    )
    plt$text(bn, size=9) / wrap_plots(plots) + plot_layout(heights=c(1,20))
}

scale_x_biexp = function() scale_x_continuous(trans=)
scale_y_biexp = function() scale_y_continuous(trans=)

fcs = list.files(pattern="\\.fcs$", recursive=TRUE, full.names=TRUE)

pdf("Rplots.pdf", 16, 10)
for (f in fcs)
    print(plot_one(f))
dev.off()

# do standard gating and subsetting here

# compensation

# bi-exponential transform

# clustering
