library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')

colors = c("no change"="#cccccc40", "down"="#006d2cd0", "up"="#a50f15d0", "different"="#045a8dd0")

stat1_cin_cor = function(go, hm) {
    ifn = bind_rows(go$stat1$wt_ifn2_over_dmso, hm$stat1$wt_ifn2_over_dmso)
    rev48 = bind_rows(go$stat1$wt_rev48_over_dmso, hm$stat1$wt_rev48_over_dmso)
    stat48 = bind_rows(go$stat1$rev48_stat1_over_wt, hm$stat1$rev48_stat1_over_wt)
    aneup = bind_rows(go$aneup, hm$aneup)

#    rev48 = hm$stat1$wt_rev48_over_dmso
#    stat48 = hm$stat1$rev48_stat1_over_wt
#    aneup = hm$aneup

    plot_one = function(df) {
        merged = inner_join(df, aneup, by="label") %>%
            mutate(type = case_when(
                abs(statistic.x) < 1.5 | abs(statistic.y) < 1.5 ~ "no change",
                statistic.x > 0 & statistic.y > 0 ~ "up",
                statistic.x < 0 & statistic.y < 0 ~ "down",
                TRUE ~ "different"
            )) %>%
            mutate(type2 = ifelse(type == "different", as.character(sign(statistic.x) > 0), type)) %>%
            group_by(type2) %>%
                mutate(score = scale(statistic.x)^2 + scale(statistic.y)^2,
                       lab = ifelse(rank(-score) <= 6 &
                                    (abs(statistic.x) + abs(statistic.y) > 10), label, NA)) %>%
            ungroup() %>%
            mutate(lab = ifelse(type == "no change", NA, lab))

        sums = merged %>%
            filter(type != "no change") %>%
            group_by(type2) %>%
            summarize(n = n())
        fet = broom::tidy(fisher.test(matrix(sums$n, ncol=2)))
        if (fet$estimate > 1) {
            odds = sprintf("%.0fx common enrichment", fet$estimate)
        } else {
            odds = sprintf("%.0fx difference enrichment", 1/fet$estimate)
        }

        ggplot(merged, aes(x=statistic.x, y=statistic.y, color=type)) +
            geom_point() +
            scale_color_manual(values=colors) +
            theme_classic() +
            geom_hline(yintercept=0, linetype="dashed", size=1.5, alpha=0.3) +
            geom_vline(xintercept=0, linetype="dashed", size=1.5, alpha=0.3) +
            geom_smooth(color="black") +
            ggrepel::geom_label_repel(aes(label=lab), size=3, max.iter=1e5, label.size=NA,
                min.segment.length=0, max.overlaps=Inf, segment.alpha=0.3, fill="#ffffffc0",
                label.padding=unit(0.2, "lines")) +
            coord_cartesian(clip="off") +
            labs(subtitle = sprintf("%s (p=%.1g FET)", odds, fet$p.value))
    }

    p1 = plot_one(ifn) + labs(title = "Acute inflammation",
                              x = "BT549: 2h IFN vs. DMSO (Wald st.)",
                              y = "Aneuploidy TPS cohort (Wald st.)")
    p2 = plot_one(rev48) + labs(title = "Acute CIN",
                                x = "BT549: 48h reversine vs. DMSO (Wald st.)",
                                y = "Aneuploidy TPS cohort (Wald st.)")
    p3 = plot_one(stat48) + labs(title = "Non-Stat1 CIN",
                                 x = "BT549 rev: 48h STAT1 KO vs. WT (Wald st.)",
                                 y = "Aneuploidy TPS cohort (Wald st.)")

    p1 + p2 + p3 + plot_layout(guides="collect")
}

sys$run({
    args = sys$cmd$parse(
        opt('g', 'cor_go', 'rds', '../expr_stat1/aneup_ko_cor/GO_Biological_Process_2020.rds'),
        opt('h', 'cor_hm', 'rds', '../expr_stat1/aneup_ko_cor/MSigDB_Hallmark_2020.rds'),
        opt('p', 'plotfile', 'pdf', 'Fig4-Stat1KO_cells.pdf')
    )

    go = readRDS(args$cor_go)
    hm = readRDS(args$cor_hm)

    expr_cor = stat1_cin_cor(go, hm)

    pdf(args$plotfile, 16, 5)
    print(expr_cor)
    dev.off()
})
