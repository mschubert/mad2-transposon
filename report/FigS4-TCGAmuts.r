library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)
sys = import('sys')
mf = import('../misc/mut_frac/mut_frac')
md = import('../misc/mut_frac/mut_distr')

mut_fracs = function() {
    mfrac = readRDS("../misc/mut_frac/mut_frac.rds")
    row1 = mf$plot_muts(mfrac$muts)
    row2 = mf$plot_cnas(mfrac$cnas)
    asm = (row1 / row2) + plot_layout(guides="collect") &
        theme_classic() & theme(
            strip.background = element_blank(),
            strip.text.x = element_text(size=10, face="bold", margin=margin(b=10))
        )
}

plot_distr = function(res) {
    qq = quantile(res$estimate, 0.7)
    res2 = res[res$gene %in% c("JAK1", "STAT1", "TTN", "TP53"),]
    p2 = ggplot(res, aes(x=estimate)) +
        geom_density() +
        geom_point(data=res2, y=0) +
        annotate("rect", xmin=qq, xmax=Inf, ymin=-Inf, ymax=Inf, fill="red", alpha=0.05) +
        geom_vline(xintercept=qq, linetype="dashed", color="grey") +
        ggrepel::geom_text_repel(data=res2, y=0, aes(label=gene), nudge_y=0.1, point.padding=2) +
        annotate("text", x=qq, y=1.85, hjust=-0.3, label="Top 30%") +
        labs(x = "Odds difference",
             y = "Density") +
        theme_minimal()
}

plot_bars = function(fet_hm) {
    hms = fet_hm %>%
        mutate(rank = dense_rank(sign(estimate) * p.value)) %>% arrange(rank) %>%
        filter(label %in% c(head(label, 11), tail(label, 9))) %>%
        mutate(label = sprintf("%s (%i)", label, size_used),
               label = forcats::fct_reorder(label, rank, .desc=TRUE)) %>%
        mutate(type = case_when(
            grepl("Spindle", label) ~ "CIN causing",
            grepl("Interf|Inflamm|Rej|TNF|STAT", label) ~ "Inflammatory",
            TRUE ~ "Other"
        )) %>%
        mutate(sign = ifelse(adj.p < 0.15, "< 0.15", "n.s."))
    ggplot(hms, aes(x=estimate, y=label, fill=type)) +
        geom_col(aes(alpha=sign)) +
        geom_vline(xintercept=1) +
        geom_text(aes(label=sprintf("  p=%.2g  ", p.value)), x=0, size=3,
            hjust=ifelse(hms$estimate > 1, 1, 0), color="grey") +
        scale_alpha_manual(values=c("< 0.15"=1, "n.s."=0.4), name="FDR") +
        scale_fill_manual(name="Gene set type",
            values=c("CIN causing"="#e41a1c", "Inflammatory"="#984ea3", "Other"="grey")) +
        scale_x_log10(breaks=c(0.5,1,2,3,5)) +
        labs(x = "Enrichment in top 30% of genes (odds ratio)") +
        theme_minimal() +
        theme(axis.title.y = element_blank())
}

sys$run({
    dens = readRDS("../misc/mut_frac/mut_distr.rds")

    p1 = fracs = mut_fracs()
    p2 = plot_distr(dens$res)
    p3 = plot_bars(dens$fet_hm)

    btm = (p2 | p3) + plot_layout(widths=c(1,1))
    asm = (wrap_plots(p1) / btm) + plot_layout(heights=c(2,1)) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("FigS4-TCGAmuts.pdf", 12, 12)
    print(asm)
    dev.off()
})
