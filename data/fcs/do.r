library(dplyr)
library(patchwork)
library(flowCore)
plt = import('plot')

ggfacs = function(df, aes, trans="identity") {
    ggplot(df, aes) +
        geom_bin2d(bins = 70) +
        scale_fill_continuous(type="viridis", trans=trans) +
        geom_density_2d(aes(color=cl), bins=7, size=0.5, linetype="dashed") +
        scale_color_brewer(palette="Set1") +
        theme_bw() +
        xlim(0, 2.5e5) +
        ylim(0, 2.5e5) +
        ggtitle(sprintf("%s vs. %s", all.vars(aes[[1]]), all.vars(aes[[2]])))
}

plot_one = function(fname, cluster=TRUE) {
    ff = read.FCS(fname)
    bn = tools::file_path_sans_ext(basename(fname))
    bn = sub("Specimen_002_", "", bn, fixed=TRUE)

    meta = ff@parameters@data
    fields = c(na.omit(meta$desc))
    df = as_tibble(ff@exprs) %>%
        dplyr::filter(`FSC-H`+`SSC-H` > 5e4, `FSC-H` > 2.5e4)
    colnames(df) = ifelse(is.na(meta$desc), meta$name, meta$desc)
#    comp = ff@description$SPILL
#    name = ff@description$`TUBE NAME`

    if (cluster) {
        mdf = df[sample(seq_len(nrow(df)), 500), fields] %>%
            data.matrix() %>% pmax(1) %>% log10()
        mcl = mclust::Mclust(mdf, G=1:5, modelNames=c("VVV"))
        pd = df[colnames(mcl$data)] %>% data.matrix() %>% pmax(1) %>% log10()
        df$cl = factor(predict(mcl, newdata=as.data.frame(pd))$classification)
    } else {
        df$cl = NA
    }

    plots = list(
        ggfacs(df, aes(x=`FSC-H`, y=`SSC-H`)),
        ggfacs(df, aes(x=`CD45`, y=`FSC-A`)) + scale_x_log10(),
        ggfacs(df, aes(x=`MAC1/GR1`, y=`SSC-A`)) + scale_x_log10(),
        ggfacs(df, aes(x=`FSC-A`, y=`B220`)) + scale_y_log10(),
        ggfacs(df, aes(x=`CD3`, y=`CD19`)) + scale_x_log10() + scale_y_log10(),
        ggfacs(df, aes(x=`SCA-1`, y=`CKIT`)) + scale_x_log10() + scale_y_log10()
    )
    plt$text(bn, size=9) / wrap_plots(plots) + plot_layout(heights=c(1,20))
}

#scale_x_biexp = function() scale_x_continuous(trans=)
#scale_y_biexp = function() scale_y_continuous(trans=)

fcs = list.files(pattern="\\.fcs$", recursive=TRUE, full.names=TRUE)

pdf("Rplots.pdf", 16, 10)
for (f in fcs)
    print(plot_one(f))
dev.off()


#fname = "./FCS files - part 2/Specimen_002_180 20_017.fcs"
#plot_one(fname)



# do standard gating and subsetting here

# compensation

# bi-exponential transform

# clustering
