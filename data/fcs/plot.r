library(patchwork)
library(dplyr)
sys = import('sys')
plt = import('plot')

#' Compute cluster centers from FACS data
#'
#' @param gated
#' @return
cluster_centers = function(gated) {
    dens_max = function(x) {
        tryCatch({
            d = density(log10(x[x>=1]), adjust=1.5) # todo: use the transformer here (or use ccs right away)
            10^(d$x[which.max(d$y)])
        }, error = function(e) NA)
    }
    gated %>%
        group_by(cl) %>%
        summarize_all(dens_max)
}

#' Plot one FACS panel
#'
#' @param df      A flowFrame-like tibble [events x markers]
#' @param meta    A tibble with fields: name [scatter, color], desc [marker]
#' @param gates
#' @param ccs     A tibble of cluster centers: cl [cluster], markers
#' @param trans   Scale transformation for axes
#' @param ctrans  Transformation for the color space (default: identity)
#' @return        A ggplot2 object of a panel
ggfacs = function(df, meta, ccs, aes, trans, gates=list(), ctrans="identity") {
    meta$desc = ifelse(is.na(meta$desc), meta$name, meta$desc)
    vs = sapply(aes, all.vars) # variable names
    cols = meta$name[match(vs, meta$desc)]
    trs = setNames(trans[cols], c("x", "y"))
    trs[sapply(trs, is.null)] = "identity"
    lims = lapply(cols, function(c) quantile(df[[c]], c(0.02, 0.98)))

    colnames(df)[match(meta$name, colnames(df))] = meta$desc
    colnames(ccs)[match(meta$name, colnames(ccs))] = meta$desc

    if (length(gates) != 0)
        gates = list(geom_polygon(data=as.data.frame(gates@boundaries), color="red", fill=NA, size=1))

    p = ggplot(df, aes) +
        geom_bin2d(bins = 70) +
        scale_fill_continuous(type="viridis", trans=ctrans) +
        geom_density_2d(aes(color=cl), bins=7, size=0.5, linetype="dashed") +
        scale_color_brewer(palette="Set1") +
        geom_text(data=ccs, aes(label="X", color=cl), size=7) +
        theme_bw() +
        scale_x_continuous(trans=trs$x, limits=lims$x) +
        scale_y_continuous(trans=trs$y, limits=lims$y) +
        gates +
        coord_cartesian(expand=FALSE) +
        labs(title = sprintf("%s vs. %s", vs[1], vs[2]),
             x = sprintf("%s [%s]", vs[1], cols[1]),
             y = sprintf("%s [%s]", vs[2], cols[2]))
    plt$try(p)
}

#' Assemble FACS panels
#'
#' @param df      A flowFrame-like tibble [events x markers]
#' @param meta    A tibble with fields: name [scatter, color], desc [marker]
#' @param gates
#' @param title
assemble = function(df, meta, trans, gates, title) {
    gated = df %>% filter(debris_gate)
    ccs = cluster_centers(gated)

    plots = list(
        ggfacs(df, meta, ccs, aes(x=`FSC-H`, y=`SSC-H`), trans, ctrans="log", gates=gates),
        ggfacs(gated, meta, ccs, aes(x=`CD45`, y=`FSC-A`), trans),
        ggfacs(gated, meta, ccs, aes(x=`MAC1/GR1`, y=`SSC-A`), trans),
        ggfacs(gated, meta, ccs, aes(x=`FSC-A`, y=`B220`), trans),
        ggfacs(gated, meta, ccs, aes(x=`CD3`, y=`CD19`), trans),
        ggfacs(gated, meta, ccs, aes(x=`SCA-1`, y=`CKIT`), trans)
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

    (plt$text(title, size=7) | pp) / wrap_plots(plots) + plot_layout(heights=c(1,20))
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'FCS files - part 1.rds'),
        opt('p', 'plotfile', 'pdf', 'FCS files - part 1.pdf')
    )

    dset = readRDS(args$infile)

    pdf(args$plotfile, 16, 10)
    for (i in seq_along(dset$res)) {
        title = names(dset$res)[i]
        message(title)
        print(with(dset, assemble(res[[i]], meta, trans, gates, title)))
    }
    dev.off()
})
