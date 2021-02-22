library(dplyr)
library(patchwork)
library(flowCore)
library(mclust)
plt = import('plot')

ggfacs = function(df, meta, ccs, aes, ctrans="identity", gate=NULL) {
    meta$desc = ifelse(is.na(meta$desc), meta$name, meta$desc)
    vs = sapply(aes, all.vars)
    scs = grepl("[FS]SC-[AH]|Time", vs) + 1 # x,y axes: 1=log, 2=linear
    logs = c("log", "identity")[scs]
    brk = list(c(10,100,1e3,1e4,1e5), 0:10 * 2.5e4)[scs]
#    xlims = quantile(df[[vs[1]]], c(0.01, 0.99))
#    ylims = quantile(df[[vs[2]]], c(0.01, 0.99))
    xlims = ylims = c(5, 2e5)
    if (is.null(gate)) {
        gates = list()
    } else {
        gates = list(geom_polygon(data=gate, color="red", fill=NA, size=1))
    }
    ggplot(df, aes) +
        geom_bin2d(bins = 70) +
        scale_fill_continuous(type="viridis", trans=ctrans) +
        geom_density_2d(aes(color=cl), bins=7, size=0.5, linetype="dashed") +
        scale_color_brewer(palette="Set1") +
        geom_text(data=ccs, aes(label="X", color=cl), size=7) +
        theme_bw() +
        scale_x_continuous(limits=range(brk[[1]]), breaks=brk[[1]], trans=logs[1]) +
        scale_y_continuous(limits=range(brk[[2]]), breaks=brk[[2]], trans=logs[2]) +
        gates +
        coord_cartesian(expand=FALSE) +
        labs(title = sprintf("%s vs. %s", vs[1], vs[2]),
             x = sprintf("%s [%s]", vs[1], meta$name[meta$desc==vs[1]]),
             y = sprintf("%s [%s]", vs[2], meta$name[meta$desc==vs[2]]))
}

log_mat = function(df, meta) { #, n_max=1.5e4) {
    fields = c(na.omit(meta$desc))
#    cl_df_idx = sample(seq_len(nrow(df)), min(nrow(df), n_max))
    mdf = df[,fields] %>%
        data.matrix() %>% pmax(1) %>% log10()
}

refine_mclust = function(df, meta, G=1:4, n_pts=500) {
    df$cl = as.character(df$cl)
    for (cl in setdiff(df$cl, NA)) { # filter by size?
        cur_cl = !is.na(df$cl) & df$cl == cl
        lmat = log_mat(df[cur_cl,], meta)
        cur_subs = sample(seq_len(nrow(lmat)), min(nrow(lmat), n_pts)) # needs 1000, G=1:5 if hdbscan finds nothing
        res = mclust::Mclust(lmat[cur_subs,], G=G, modelNames="VVV") # 1:4 only if HDB didn't find a cluster?
        df$cl[cur_cl] = paste0(cl, predict(res, newdata=lmat)$classification)
    }
    df
}

plot_one = function(fname, cluster=TRUE) {
    message(fname)
    ff = read.FCS(fname)
    bn = tools::file_path_sans_ext(basename(fname))
    bn = sub("Specimen_002_", "", bn, fixed=TRUE)

    gates = yaml::read_yaml("fsc-ssc-gates.yaml")
    fsc_ssc = as_tibble(gates$common$debris) * 1e3
    if (basename(fname) %in% names(gates$sample))
        fsc_ssc = as_tibble(gates$sample[[basename(fname)]]$debris) * 1e3
    debris = do.call(polygonGate, c(fsc_ssc, list(filterId="debris")))

    meta = ff@parameters@data
    fields = c(na.omit(meta$desc))
    ungated = as_tibble(ff@exprs)
    db = flowCore::filter(ff, debris)
    df = as_tibble(ff@exprs)[db@subSet,] %>%
        dplyr::filter(`SSC-A` < 2.5e5) # todo: is this caught by singlets [FSC-H FSC-A] gate?
    df = df[rowSums(df < 0) == 0,]
    colnames(df) = ifelse(is.na(meta$desc), meta$name, meta$desc)
    colnames(ungated) = ifelse(is.na(meta$desc), meta$name, meta$desc)
#    comp = ff@description$SPILL
#    name = ff@description$`TUBE NAME`

    cl_df_idx = sample(seq_len(nrow(df)), min(nrow(df), 1.5e4))
    mdf = df[cl_df_idx, fields] %>%
        data.matrix() %>% pmax(1) %>% log10()
    mcl = dbscan::hdbscan(mdf, minPts=200)
    G = 1:5
    mclust_pts = 1500 # e.g. 422 not split ckit +/- with 1000 pts
    if (length(unique(mcl$cluster)) > 1) {
        mcl$cluster[mcl$cluster == 0] = NA
        G = 1:4
        mclust_pts = 500
    }
    df$cl = NA
    df$cl[cl_df_idx] = mcl$cluster
    df = refine_mclust(df, meta, G, mclust_pts)
    df$cl = factor(df$cl)
    ungated = left_join(ungated, df) # add cluster information
#    levels(df$cl) = seq_along(levels(df$cl))

    dens_max = function(x) {
        d = density(log10(x[x>=1]), adjust=1.5)
        10^(d$x[which.max(d$y)])
    }
    ccs = df %>%
        group_by(cl) %>%
        summarize_all(dens_max)

    plots = list(
        ggfacs(ungated, meta, ccs, aes(x=`FSC-H`, y=`SSC-H`), ctrans="log", gate=fsc_ssc),
        ggfacs(df, meta, ccs, aes(x=`CD45`, y=`FSC-A`)),
        ggfacs(df, meta, ccs, aes(x=`MAC1/GR1`, y=`SSC-A`)),
        ggfacs(df, meta, ccs, aes(x=`FSC-A`, y=`B220`)),
        ggfacs(df, meta, ccs, aes(x=`CD3`, y=`CD19`)),
        ggfacs(df, meta, ccs, aes(x=`SCA-1`, y=`CKIT`))
    )

    not_na = sum(!is.na(df$cl))
    fdf = group_by(df, cl) %>%
        summarize(pro = n() / not_na) %>%
        mutate(lab = sprintf("(%i) %.0f%%", cl, pro * 100))
    pp = ggplot(fdf) +
        geom_point(aes(x=cl, size=pro, color=factor(cl)), y=1) +
        scale_color_brewer(palette="Set1") +
        geom_text(aes(x=cl, label=lab), y=0.5) +
        scale_size_area(guide=FALSE) +
        theme_void() +
        theme(legend.position = "none")

    list(
        ccs = ccs,
        plot = (plt$text(bn, size=7) | pp) / wrap_plots(plots) + plot_layout(heights=c(1,20))
    )
}

#scale_x_biexp = function() scale_x_continuous(trans=)
#scale_y_biexp = function() scale_y_continuous(trans=)

dir = "FCS files - part 1"
fcs = list.files(dir, pattern="\\.fcs$", recursive=TRUE, full.names=TRUE)

#fcs = c(
#    "FCS files - part 1/Specimen_002_401 21_019.fcs",
#    "FCS files - part 1/Specimen_002_443 20_046.fcs"
#)

res = sapply(fcs, function(x) try(plot_one(x)), simplify=FALSE)
errs = sapply(res, class) == "try-error"
if (any(errs)) {
    warning("Samples failed: ", paste(names(res)[errs], collapse=", "))
    res = res[!errs]
}

ccs = lapply(res, function(x) x$ccs)
plots = lapply(res, function(x) x$plot)

pdf("Rplots.pdf", 16, 10)
for (p in plots)
    print(plt$try(p))
dev.off()

saveRDS(ccs, file="ccs.rds")


#fname = "./FCS files - part 2/Specimen_002_180 20_017.fcs"
#fname = "FCS files - part 1/Specimen_002_404 20_022.fcs"
#fname = "FCS files - part 1/Specimen_002_404 20_022.fcs"
#plot_one(fname)
fname = "FCS files - part 1/Specimen_002_446 20_048.fcs"
fname = "FCS files - part 1/Specimen_002_443 20_046.fcs"
fname = "FCS files - part 1/Specimen_002_422 20_033.fcs" # 1500 pts split ckit+/-?

# do standard gating and subsetting here

# compensation

# bi-exponential transform

# clustering
