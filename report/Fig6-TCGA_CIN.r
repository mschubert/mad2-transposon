library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(patchwork)
library(ggpmisc)
sys = import('sys')

survplot = function(dset) {
    p53wt = dset %>% filter(p53_mut == 0) #, `Myc Targets V1`>0)
    p53mut = dset %>% filter(p53_mut != 0)

    m1 = coxph(Surv(OS_years, vital_status) ~ age_at_diagnosis + purity + iclass, data=p53wt)
    m1p = broom::tidy(m1) %>% filter(term == "iclassCIN_stat1ko") %>% pull(p.value)
    m2 = coxph(Surv(OS_years, vital_status) ~ age_at_diagnosis + purity + iclass, data=p53mut)
    m2p = broom::tidy(m2) %>% filter(term == "iclassCIN_stat1ko") %>% pull(p.value)

    pal = c("#ababab", "blue", "#ad07e3")
    lab = c("No CIN signature", "Acute CIN signature", "STAT1ko CIN signature")

    fit1 = survfit(Surv(OS_years, vital_status) ~ iclass, data=p53wt)
    ps1 = ggsurvplot(fit1, data=p53wt, xlim=c(0,10), break.time.by=2.5, palette=pal, legend.labs=lab)$plot +
        ylim(c(0.25,1)) +
        labs(x="Overall survival (years)", title="Acute CIN response vs. STAT1 ko", subtitle="p53 wt") +
        annotate("text_npc", npcx=0.1, npcy=0.1,
                 label=sprintf("CIN STAT1ko vs. CIN p=%.2g", m1p))

    fit2 = survfit(Surv(OS_years, vital_status) ~ iclass, data=p53mut)
    ps2 = ggsurvplot(fit2, data=p53mut, xlim=c(0,10), break.time.by=2.5, palette=pal, legend.labs=lab)$plot +
        ylim(c(0.25,1)) +
        labs(x="Overall survival (years)", subtitle="p53 mut") +
        annotate("text_npc", npcx=0.1, npcy=0.1,
                 label=sprintf("CIN STAT1ko vs. CIN p=%.2g", m2p))

    ps1 + ps2 + plot_layout(guides="collect") & theme(legend.direction = "vertical")
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
    anov = list(MYC = aov_one(MYC ~ aneup_log2seg +type+p53_mut + myc_copy + rev48_stat1_over_wt),
                MycV1 = aov_one(`Myc Targets V1` ~ aneup_log2seg +type+p53_mut + MYC + myc_copy + rev48_stat1_over_wt),
                aneup_log2seg = aov_one(aneup_log2seg ~ `Myc Targets V1`+type + p53_mut + MYC + rev48_stat1_over_wt)) %>%
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
              panel.background = element_blank())
}

purplot = function(dset) {
    dset %>%
        mutate(p53_mut = as.character(p53_mut), #TODO: check if wt/mut ordered right
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
        ), iclass=relevel(factor(iclass), "noCIN")) %>%
        mutate(mclass = case_when(
            `Myc Targets V1`< 0 ~ "nomyc",
            rev48_stat1_over_wt>0 ~ "stat1ko_myc",
            rev48_stat1_over_wt<0 ~ "stat1wt_myc"
        ), mclass=relevel(factor(mclass), "nomyc"))
})

# overview BRCA cohort

# CIN surv + infiltration
