library(dplyr)
library(patchwork)
library(flowCore)
plt = import('plot')

ggfacs = function(ff, aes, mcl=NULL, trans="identity") {
    meta = ff@parameters@data
    df = as_tibble(ff@exprs) %>%
        dplyr::filter(`FSC-H`+`SSC-H` > 5e4, `FSC-H` > 2.5e4)
    colnames(df) = ifelse(is.na(meta$desc), meta$name, meta$desc)
    comp = ff@description$SPILL
    name = ff@description$`TUBE NAME`

    if (!is.null(mcl)) {
        pd = log10(data.matrix(df[colnames(mcl$data)]))
        pd[is.na(pd)] = 0
        pd[is.infinite(pd)] = 0
        df$cl = factor(predict(mcl, newdata=as.data.frame(pd))$classification)
    }

    ggplot(df, aes) +
        geom_bin2d(bins = 70) +
        scale_fill_continuous(type="viridis", trans=trans) +
        geom_density_2d(aes(color=cl), bins=7, size=0.25) +
        theme_bw() +
        xlim(0, 2.5e5) +
        ylim(0, 2.5e5) +
        ggtitle(sprintf("%s vs. %s", all.vars(aes[[1]]), all.vars(aes[[2]])))
}

plot_one = function(fname, mcl=NULL) {
    fc = read.FCS(fname)
    bn = tools::file_path_sans_ext(basename(fname))
    bn = sub("Specimen_002_", "", bn, fixed=TRUE)
    plots = list(
        ggfacs(fc, aes(x=`FSC-H`, y=`SSC-H`), mcl),
        ggfacs(fc, aes(x=`CD45`, y=`FSC-A`), mcl) + scale_x_log10(),
        ggfacs(fc, aes(x=`MAC1/GR1`, y=`SSC-A`), mcl) + scale_x_log10(),
        ggfacs(fc, aes(x=`FSC-A`, y=`B220`), mcl) + scale_y_log10(),
        ggfacs(fc, aes(x=`CD3`, y=`FSC-A`), mcl) + scale_x_log10(),
        ggfacs(fc, aes(x=`SCA-1`, y=`CKIT`), mcl) + scale_x_log10() + scale_y_log10()
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
fields = c("CKIT","CD3","B220","MAC1/GR1","CD19","CD45","SCA-1")
xorig = df[sample(seq_len(nrow(df)), 500),]
#xorig = df[,fields]
x = log10(data.matrix(xorig[,fields]))
x[is.na(x)] = 0
x[is.infinite(x)] = 0
#x = apply(x, 2, scale)


#umap2 = uwot::umap(x, n_neighbors=500, n_components=3, scale=TRUE)
#cl = igraph::cluster_louvain(scran::buildSNNGraph(t(umap2), k=10))$membership
#table(cl)

library(mclust)
res = Mclust(x, modelNames=c("VVV"))
#library(dbscan)
#cl = dbscan(x, eps=1, minPts=50)




plot_one(fname, res)
# how about umap? (+cluster, mean+sd in original space, annotate)




# do standard gating and subsetting here

# compensation

# bi-exponential transform

# clustering
