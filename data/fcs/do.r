library(dplyr)
library(patchwork)
library(flowCore)
plt = import('plot')

ggfacs = function(ff, aes, ccs=data.frame(), trans="identity") {
    meta = ff@parameters@data
    df = as_tibble(ff@exprs) %>%
        dplyr::filter(`FSC-H`+`SSC-H` > 5e4, `FSC-H` > 2.5e4)
    colnames(df) = ifelse(is.na(meta$desc), meta$name, meta$desc)
    comp = ff@description$SPILL
    name = ff@description$`TUBE NAME`

    for (cn in setdiff(colnames(df), colnames(ccs)))
        ccs[[cn]] = NA_real_

    ggplot(df, aes) +
        geom_bin2d(bins = 70) +
        scale_fill_continuous(type="viridis", trans=trans) +
        geom_text(data=ccs, aes(label="X", color=cl), size=10) +
        theme_bw() +
        xlim(0, 2.5e5) +
        ylim(0, 2.5e5) +
        ggtitle(sprintf("%s vs. %s", all.vars(aes[[1]]), all.vars(aes[[2]])))
}

plot_one = function(fname, ccs) {
    fc = read.FCS(fname)
    bn = tools::file_path_sans_ext(basename(fname))
    bn = sub("Specimen_002_", "", bn, fixed=TRUE)
    plots = list(
        ggfacs(fc, aes(x=`FSC-H`, y=`SSC-H`), ccs),
        ggfacs(fc, aes(x=`CD45`, y=`FSC-A`), ccs) + scale_x_log10(),
        ggfacs(fc, aes(x=`MAC1/GR1`, y=`SSC-A`), ccs) + scale_x_log10(),
        ggfacs(fc, aes(x=`FSC-A`, y=`B220`), ccs) + scale_y_log10(),
        ggfacs(fc, aes(x=`CD3`, y=`FSC-A`), ccs) + scale_x_log10(),
        ggfacs(fc, aes(x=`SCA-1`, y=`CKIT`), ccs) + scale_x_log10() + scale_y_log10()
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


fname = "./FCS files - part 2/Specimen_002_180 20_017.fcs"
ff = read.FCS(fname)
    meta = ff@parameters@data
    df = as_tibble(ff@exprs) %>%
        dplyr::filter(`FSC-H`+`SSC-H` > 5e4, `FSC-H` > 2.5e4)
    colnames(df) = ifelse(is.na(meta$desc), meta$name, meta$desc)
plot_one(fname)
fields = c("FSC-A","SSC-A","CKIT","CD3","B220","MAC1/GR1","CD19","CD45","SCA-1")
xorig = df[sample(seq_len(nrow(df)), 1e4), fields]
x = data.matrix(xorig)
x[,-c(1:2)] = log10(x[,-c(1:2)])
x[is.na(x)] = 0
x[is.infinite(x)] = 0
#library(mclust)
#res = Mclust(x, G=4)

library(dbscan)
cl = hdbscan(x, minPts=100)
ccs = xorig %>%
    mutate(cl = factor(cl$cluster)) %>%
    group_by(cl) %>%
    mutate(n = n()) %>%
    summarize_all(median) %>%
    dplyr::filter(cl != "0")

# how about umap? (+cluster, mean+sd in original space, annotate)




# do standard gating and subsetting here

# compensation

# bi-exponential transform

#cc clustering
