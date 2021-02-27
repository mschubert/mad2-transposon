library(patchwork)
library(dplyr)
sys = import('sys')
plt = import('plot')

ggfacs = function(df, meta, ccs, aes, ctrans="identity", gate=NULL) {
    df = df %>% filter(CKIT >= 1) # neg due compensation, needs biexp trans instead log
#    meta$desc = ifelse(is.na(meta$desc), meta$name, meta$desc)
    vs = sapply(aes, all.vars)
    scs = grepl("[FS]SC-[AH]|Time", vs) + 1 # x,y axes: 1=log, 2=linear
    logs = list("log", "identity")[scs]
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
             x = sprintf("%s [%s]", vs[1], vs[1]), #meta$name[meta$desc==vs[1]]),
             y = sprintf("%s [%s]", vs[2], vs[2]))#, meta$name[meta$desc==vs[2]]))
}

plot_one = function(df, title) {
    message(title)
    # would be better to save this with object
    gates = yaml::read_yaml("fsc-ssc-gates.yaml")
    fsc_ssc = as_tibble(gates$common$debris) * 1e3
#    if (basename(fname) %in% names(gates$sample))
#        fsc_ssc = as_tibble(gates$sample[[basename(fname)]]$debris) * 1e3

    df$cl = factor(df$cl) #todo: in data
    gated = df %>% filter(debris_gate)

    dens_max = function(x) {
        d = density(log10(x[x>=1]), adjust=1.5)
        10^(d$x[which.max(d$y)])
    }
    ccs = gated %>%
        group_by(cl) %>%
        summarize_all(dens_max)

    plots = list(
        ggfacs(df, meta, ccs, aes(x=`FSC-H`, y=`SSC-H`), ctrans="log", gate=fsc_ssc),
        ggfacs(gated, meta, ccs, aes(x=`CD45`, y=`FSC-A`)),
        ggfacs(gated, meta, ccs, aes(x=`MAC1/GR1`, y=`SSC-A`)),
        ggfacs(gated, meta, ccs, aes(x=`FSC-A`, y=`B220`)),
        ggfacs(gated, meta, ccs, aes(x=`CD3`, y=`CD19`)),
        ggfacs(gated, meta, ccs, aes(x=`SCA-1`, y=`CKIT`))
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

    res = readRDS(args$infile)
    plots = mapply(function(...) try(plot_one(...)), res, names(res), SIMPLIFY=FALSE)

    plots = plots[sapply(plots, class) != "try-error"]

    pdf(args$plotfile, 16, 10)
    for (p in plots)
        try(print(plt$try(p)))
    dev.off()
})
