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

plot_bars = function(fet_hm) {
    hms = fet_hm %>%
        filter(label %in% c(head(label, 10), tail(label, 10))) %>%
        mutate(label = forcats::fct_reorder(label, p.value, .desc=TRUE)) %>%
        mutate(type = case_when(
            grepl("Spindle", label) ~ "CIN causing",
            grepl("Interf|Inflamm|Rej|TNF|STAT", label) ~ "Inflammatory",
            TRUE ~ "Other"
        )) %>%
        mutate(sign = ifelse(adj.p < 0.15, "< 0.15", "n.s."))
    ggplot(hms, aes(x=estimate, y=label, fill=type)) +
        geom_col(aes(alpha=sign)) +
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

    qq = quantile(dens$res$estimate, 0.7)
    p2 = md$distr_plot(dens$res) +
        annotate("rect", xmin=qq, xmax=Inf, ymin=-Inf, ymax=Inf, fill="red", alpha=0.05) +
        geom_vline(xintercept=qq, linetype="dashed") +
        annotate("text", x=qq, y=1.85, hjust=-0.3, label="Top 30%") +
        labs(x = "Odds difference",
             y = "Density") +
        theme_minimal()
    p3 = plot_bars(dens$fet_hm) + geom_vline(xintercept=1)

    btm = (p2 | p3) + plot_layout(widths=c(1,1))
    asm = (wrap_plots(p1) / btm) + plot_layout(heights=c(2,1)) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf("FigS4-TCGAmuts.pdf", 12, 12)
    print(asm)
    dev.off()
})
