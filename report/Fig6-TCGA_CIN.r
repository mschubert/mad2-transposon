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

surv_atn = function(dset) {
    # atn explanation for BRCA surv assocs (E2F, MycV1, CIN70, aneup, [etc.])
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

    survplot(dset)
})
