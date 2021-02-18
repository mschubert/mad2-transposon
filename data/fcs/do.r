library(dplyr)
library(patchwork)
library(flowCore)
library(flowViz)

ggfacs = function(ff, aes) {
    meta = ff@parameters@data
    df = as_tibble(ff@exprs)
    colnames(df) = ifelse(is.na(meta$desc), meta$name, meta$desc)
    comp = ff@description$SPILL
    name = ff@description$`TUBE NAME`

    ggplot(df, aes) +
        geom_bin2d(bins = 70) +
        scale_fill_continuous(type = "viridis") +
        theme_bw() +
        xlim(0, 2.5e5) +
        ylim(0, 2.5e5) +
        ggtitle(sprintf("%s vs. %s", all.vars(aes[[1]]), all.vars(aes[[2]])))
}

fcs = list.files(pattern="\\.fcs$", recursive=TRUE, full.names=TRUE)

fc = read.FCS("./FCS files - part 2/Specimen_002_180 20_017.fcs") # testing
plot(fc, c("FSC-A", "SSC-A"), nbin=100)
plot(fc, c("FSC-A", "PerCP-Cy5-5-A"), nbin=100) + ggtitle("test") # CD45
plot(fc, c("FSC-A", "SSC-A"), nbin=100)

plots = list(
    ggfacs(fc, aes(x=`FSC-H`, y=`SSC-H`)),
    ggfacs(fc, aes(x=`FSC-H`, y=`CD45`)) + scale_y_log10(),
    ggfacs(fc, aes(x=`MAC1/GR1`, y=`SSC-A`)) + scale_x_log10(),
    ggfacs(fc, aes(x=`FSC-A`, y=`B220`)) + scale_y_log10(),
    ggfacs(fc, aes(x=`CD3`, y=`FSC-H`)) + scale_x_log10(),
    ggfacs(fc, aes(x=`SCA-1`, y=`CKIT`)) + scale_x_log10() + scale_y_log10()
)
wrap_plots(plots)

# do standard gating and subsetting here

# compensation

# bi-exponential transform

# clustering
