library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
sys = import('sys')

isCINsig_plot = function(dset) {
    cins = tidyr::gather(dset, "measure", "value", aneup_log2seg, CIN70_Carter2006, HET70)
    stats = cins %>% group_by(measure) %>%
        summarize(res = list(lm(value ~ rev48_stat1_over_wt) %>% broom::tidy() %>%
                             filter(term == "rev48_stat1_over_wt"))) %>%
        tidyr::unnest_wider(res)

    ggplot(cins, aes(x=rev48_stat1_over_wt, y=value)) +
        geom_point(aes(fill=type, shape=factor(p53_mut)), size=2, alpha=0.6) +
        geom_smooth(method="lm", se=FALSE) +
        facet_wrap(~measure, scales="free") +
        geom_text_npc(data=stats, aes(label=sprintf("p=%.2g", p.value)),
                      npcx=0.08, npcy=0.95, size=4) +
        scale_shape_manual(values=c("0"=21, "1"=23), name="TP53 mutation")
}

sgl_plot = function(dset) {
    sgl_one = function(x, y) {
        fml = as.formula(sprintf("`%s` ~ `%s`", y, x))
        list(all = lm(fml, data=dset) %>% broom::tidy(),
             wt = lm(fml, data=dset %>% filter(p53_mut == 0)) %>% broom::tidy(),
             mut = lm(fml, data=dset %>% filter(p53_mut == 1)) %>% broom::tidy()) %>%
        bind_rows(.id="p53_status") %>%
        filter(grepl(x, term, fixed=TRUE))
    }
    sgl = expand.grid(x = c("type", "p53_mut", "myc_copy", "MYC", "rev48_stat1_over_wt"),
                      y = c("MYC", "Myc Targets V1"), stringsAsFactors=FALSE) %>%
        filter(x != y) %>%
        rowwise() %>%
            mutate(res = list(sgl_one(x, y))) %>%
        ungroup() %>%
        tidyr::unnest(res) %>%
        filter(!is.na(estimate)) %>%
        mutate(p53_status = factor(p53_status, levels=c("all", "wt", "mut")))

    ggplot(sgl, aes(x=term, y=-log10(p.value))) +
        geom_col(aes(fill=factor(sign(estimate)))) +
        facet_grid(y ~ p53_status, scales="free_x", space="free_x") +
        theme(axis.text.x = element_text(hjust=1, angle=45)) +
        labs(title="Single associations", x="Predictor", fill="Regression slope")
}

aov_plot = function(dset) {
    aov_one = function(fml) {
        do_aov = function(mod) car::Anova(mod) %>% broom::tidy() %>%
            mutate(frac_sq = sumsq / sum(rev(sumsq)[-1], na.rm=TRUE),
                   frac_sq_total = sumsq / sum(sumsq, na.rm=TRUE),
                   total_sq = sum(sumsq, na.rm=TRUE) / nobs(mod))
        list(all = lm(fml, dset) %>% do_aov(),
             wt = lm(fml, dset %>% filter(p53_mut == 0)) %>% do_aov(),
             mut = lm(fml, dset %>% filter(p53_mut == 1)) %>% do_aov()) %>%
        bind_rows(.id = "p53_status") %>%
            mutate(rel_sq = total_sq / max(total_sq, na.rm=TRUE))
    }
    anov = list(MYC = aov_one(MYC ~ type + p53_mut + myc_copy + rev48_stat1_over_wt),
                MycV1 = aov_one(`Myc Targets V1` ~ type + p53_mut + MYC + myc_copy + rev48_stat1_over_wt)) %>%
        bind_rows(.id="y") %>%
        mutate(p53_status = factor(p53_status, levels=c("all", "wt", "mut")),
               y = factor(y, levels=c("MYC", "MycV1", "aneup_log2seg")))

    ggplot(anov, aes(x=rel_sq/2, fill=term, y=frac_sq_total, width=rel_sq)) +
        geom_col() +
        coord_polar("y", start=0) +
        facet_grid(y ~ p53_status) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              panel.background = element_blank()) +
        labs(title="Variance explained (Type II ANOVA)", fill="Predictor")
}

purplot = function(dset) {
    ds = dset %>%
        mutate(p53_mut = as.character(p53_mut),
               MYCg = cut(MYC, c(-Inf,median(MYC),Inf)),
               MycV1 = cut(`Myc Targets V1`, c(-Inf,0,Inf)))

    p1 = ggplot(ds, aes(x=iclass, fill=iclass, y=1-purity)) +
        geom_boxplot(outlier.shape=NA) + facet_grid(MYCg ~ p53_mut) +
        ggtitle("MYC")
    p2 = ggplot(ds, aes(x=iclass, fill=iclass, y=1-purity)) +
        geom_boxplot(outlier.shape=NA) + facet_grid(MycV1 ~ p53_mut) +
        ggtitle("Myc Targets V1")
    p1 + p2 + plot_layout(guides="collect")
}

sys$run({
    args = sys$cmd$parse(
        opt('b', 'brca', 'rds', '../tcga_myc/dset.rds'),
        opt('p', 'plotfile', 'pdf', 'FigS5-CINsupp.pdf')
    )

    brca = readRDS(args$brca)
    brca$meta$vital_status = as.integer(brca$meta$vital_status) - 1
    dset = cbind(brca$meta, as.data.frame(brca$dmat)) %>%
        mutate(OS_years = OS_time / 365,
               vital_status = ifelse(OS_years > 10, 0, vital_status),
               OS_years = pmin(OS_years, 10)) %>%
        mutate(iclass = case_when(
            wt_rev48_over_dmso<0 & rev48_stat1_over_wt>0 ~ "CIN_stat1ko",
            wt_rev48_over_dmso>0 ~ "CIN",
            wt_rev48_over_dmso<0 & rev48_stat1_over_wt<0 ~ "noCIN"
        ), iclass=relevel(factor(iclass), "noCIN"))

    p1 = isCINsig_plot(dset)
    p2 = purplot(dset)
    right = sgl_plot(dset) + aov_plot(dset) + plot_layout(nrow=2, heights=c(2,3))
    left = p1 + p2 + plot_layout(heights=c(1,2), ncol=1)
    p = wrap_plots(left) + wrap_plots(right) + plot_layout(widths=c(1.6,1)) +
        plot_annotation(tag_levels='a') &
        theme(plot.tag = element_text(size=18, face="bold"))

    pdf(args$plotfile, 16, 9)
    print(p)
    dev.off()
})
