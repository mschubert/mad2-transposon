library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(patchwork)
library(ggpmisc)
sys = import('sys')
tcga = import('data/tcga')

#todo: add what is now Fig4 BRCA sig compare

#todo: if instead of stat1ko, put in aneup/CIN70/E2F? @supp
# + naive assocs with those @supp

pancan_myc_stat = function() {
    pur = tcga$purity() %>% select(Sample, cohort, purity=estimate)
    cohorts = unique(pur$cohort)
    myc = lapply(cohorts, function(c) tcga$gsva(c, "MSigDB_Hallmark_2020")["Myc Targets V1",]) %>%
        do.call(c, .) %>% stack() %>% as_tibble() %>% select(Sample=ind, MycV1=values)
    stat = lapply(cohorts, function(c) tcga$gsva(c, "DoRothEA")["STAT1 (a)",]) %>%
        do.call(c, .) %>% stack() %>% as_tibble() %>% select(Sample=ind, STAT1=values)

    ds = tcga$aneuploidy() %>% select(Sample, aneup=aneup_log2seg) %>%
        filter(substr(Sample, 14, 16) == "01A") %>%
        inner_join(pur) %>%
        inner_join(myc) %>%
        inner_join(stat) %>%
        mutate(aneup = aneup / purity)

    x= ds %>% group_by(cohort) %>%
        summarize(res = list(broom::tidy(lm(STAT1 ~ purity + aneup)))) %>%
        tidyr::unnest(res)
    y= ds %>% group_by(cohort) %>%
        summarize(res = list(broom::tidy(lm(MycV1 ~ purity + aneup)))) %>%
        tidyr::unnest(res)

    ds2 = inner_join(
        x %>% filter(term == "aneup") %>% select(cohort, STAT1=statistic),
        y %>% filter(term == "aneup") %>% select(cohort, MycV1=statistic)
    )
    ggplot(ds2, aes(x=STAT1, y=MycV1, color=cohort)) + geom_point() +
        geom_text(aes(label=cohort))
}

calc_surv = function(dset) {
    calc_one = function(field) {
        fs = rlang::sym(field)
        dset2 = dset %>%
            mutate(field = case_when(
            wt_rev48_over_dmso<0 & !! fs>0 ~ "CIN_field",
            wt_rev48_over_dmso>0 ~ "CIN",
            wt_rev48_over_dmso<0 & !! fs<0 ~ "noCIN"
        ), field=relevel(factor(field), "noCIN"))
        fml = Surv(OS_years, vital_status) ~ age_at_diagnosis + purity + field
        list(#all = coxph(fml, data=dset2),
             wt = coxph(fml, data=dset2 %>% filter(p53_mut == 0)),
             mut = coxph(fml, data=dset2 %>% filter(p53_mut == 1))) %>%
            lapply(broom::tidy) %>% bind_rows(.id="p53_status")
    }

    fields = c("rev48_stat1_over_wt", "CIN70_Carter2006", "HET70",
               "Myc Targets V1", "E2F Targets")
    sapply(fields, calc_one, simplify=FALSE) %>% bind_rows(.id="assoc") %>%
        filter(term == "fieldCIN_field") %>% select(-term)
}

survplot = function(dset) {
    p53wt = dset %>% filter(p53_mut == 0)
    p53mut = dset %>% filter(p53_mut != 0)

    iclass_cmp = paste0("iclass", rev(levels(dset$iclass))[1])
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

    other = calc_surv(dset)
    ggplot(other, aes(x=assoc, fill=p53_status, y=-log10(p.value), alpha=p.value<0.05)) +
        geom_col(position="dodge") +
        geom_hline(yintercept=-log10(0.05), linetype="dashed") +
        geom_text(aes(label=sprintf("  %.2g  ", p.value)), position=position_dodge(width=1),
                  hjust=ifelse(other$p.value<0.3, 1, 0)) +
        scale_alpha_manual(values=c("TRUE"=0.9, "FALSE"=0.4)) +
        coord_flip()

    ps1 + ps2 + plot_layout(guides="collect") & theme(legend.direction = "vertical")
}

sys$run({
    args = sys$cmd$parse(
        opt('b', 'brca', 'rds', '../tcga_myc/dset.rds'),
        opt('p', 'plotfile', 'pdf', 'Fig6-TCGA_CIN.pdf')
    )

    scde = readRDS("../data/scRNA_cancer/dset.rds")

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
