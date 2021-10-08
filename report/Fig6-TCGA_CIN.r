library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(patchwork)
library(ggpmisc)
sys = import('sys')

#todo: if instead of stat1ko, put in aneup/CIN70/E2F? @supp
# + naive assocs with those @supp
survplot = function(dset, iclass_cmp="iclassCIN_stat1ko") {
    p53wt = dset %>% filter(p53_mut == 0)
    p53mut = dset %>% filter(p53_mut != 0)

    m1 = coxph(Surv(OS_years, vital_status) ~ age_at_diagnosis + purity + iclass, data=p53wt)
    m1p = broom::tidy(m1) %>% filter(term == iclass_cmp) %>% pull(p.value)
    m2 = coxph(Surv(OS_years, vital_status) ~ age_at_diagnosis + purity + iclass, data=p53mut)
    m2p = broom::tidy(m2) %>% filter(term == iclass_cmp) %>% pull(p.value)

    pal = c("#ababab", "blue", "#ad07e3")
    lab = c("No CIN signature", "Acute CIN signature", "STAT1ko CIN signature")

    fit1 = survfit(Surv(OS_years, vital_status) ~ iclass, data=p53wt)
    ps1 = ggsurvplot(fit1, data=p53wt, xlim=c(0,10), break.time.by=2.5, palette=pal, legend.labs=lab)$plot +
        ylim(c(0.25,1)) +
        labs(x = "Overall survival (years)",
             title = "Acute CIN response vs. STAT1 ko",
             subtitle = sprintf("p53 wt (n=%i)", sum(fit1$n))) +
        annotate("text_npc", npcx=0.1, npcy=0.1,
                 label=sprintf("CIN STAT1ko vs. no CIN p=%.2g\nCIN70 n.s.\nMyc Targets V1 n.s.\nE2F Targets n.s.", m1p))

    fit2 = survfit(Surv(OS_years, vital_status) ~ iclass, data=p53mut)
    ps2 = ggsurvplot(fit2, data=p53mut, xlim=c(0,10), break.time.by=2.5, palette=pal, legend.labs=lab)$plot +
        ylim(c(0.25,1)) +
        labs(x = "Overall survival (years)",
             subtitle = sprintf("p53 mut (n=%i)", sum(fit2$n))) +
        annotate("text_npc", npcx=0.1, npcy=0.1,
                 label=sprintf("CIN STAT1ko vs. no CIN p=%.2g\nCIN70 n.s.\nMyc Targets V1 n.s.\nE2F Targets n.s.", m2p))

    ps1 + ps2 + plot_layout(guides="collect") & theme(legend.direction = "vertical")
}

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
        theme(axis.text.x = element_text(hjust=1, angle=30)) +
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

    # todo: ideally, add survival
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
    dset %>%
        mutate(p53_mut = as.character(p53_mut), #todo: add myc expr
               MycV1 = cut(`Myc Targets V1`, c(-Inf,0,Inf))) %>%
        ggplot(aes(x=iclass, fill=iclass, y=1-purity)) +
            geom_boxplot(outlier.shape=NA) + facet_grid(MycV1 ~ p53_mut)
}

sys$run({
    args = sys$cmd$parse(
        opt('b', 'brca', 'rds', '../tcga_myc/dset.rds'),
        opt('p', 'plotfile', 'pdf', 'Fig6-TCGA_CIN.pdf')
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

    sgl_plot(dset) + aov_plot(dset) + plot_layout(nrow=2, heights=c(2,3))
})
