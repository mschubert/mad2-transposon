# data frame with:
# * myc copy number ,expression level, target expressioin
# * immune cells
# * subtype
# * survival
library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
gu = import('../tcga_stat1/genenet_util')

lm_plot = function(dset, aes, covar) {
    v = sapply(aes, all.vars) %>%
        match(colnames(dset))
    m1 = broom::tidy(lm(dset[[v[[2]]]] ~ dset[[v[[1]]]], data=dset)) %>%
        filter(term == "dset[[v[[1]]]]")
    m2 = broom::tidy(lm(dset[[v[[2]]]] ~ dset[[covar]] + dset[[v[[1]]]], data=dset)) %>%
        filter(term == "dset[[v[[1]]]]")
    ggplot(dset, aes) +
        geom_point() +
        geom_smooth(method="lm") +
        labs(subtitle = sprintf("%.2f p=%.2g (%.2f p=%.2g with %s)",
                                m1$estimate, m1$p.value, m2$estimate, m2$p.value, covar))
}

check_surv = function(var) {
    var = dset[["Myc Targets V2"]]
    survival::coxph(survival::Surv(OS_time, as.integer(vital_status)-1) ~
                    year_of_birth + sex + purity + var, data=dset) %>%
        broom::tidy()
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'rds', 'dset.rds'),
        opt('o', 'outfile', 'rds', 'brca.rds'),
        opt('p', 'plotfile', 'pdf', 'brca.pdf')
    )

    ds = readRDS(args$infile)
#    ds$meta$OS_time = ifelse(ds$meta$OS_time > 365*5, NA, ds$meta$OS_time)
    ds$meta$vital_status = as.integer(ds$meta$vital_status) - 1
    dset = cbind(ds$meta, as.data.frame(ds$dmat))

    dens = dset %>%
        select(Sample, purity, aneup_abs, aneup_log2seg, `Interferon Gamma Response`,
               `Myc Targets V1`, CIN70_Carter2006,
               wt_rev24_over_dmso, wt_rev48_over_dmso, rev24_stat1_over_wt, rev48_stat1_over_wt) %>%
        tidyr::gather("field", "value", -Sample)

    p1 = ggplot(dens, aes(x=value)) +
        geom_density() +
        facet_wrap(~ field, scales="free")

    x = na.omit(data.matrix(dset[-1]))
    xsub = c("purity", "aneup_log2seg", "Myc Targets V1", "Interferon Gamma Response",
             "wt_rev24_over_dmso", "rev24_stat1_over_wt", #"ifn_loss",
             "myc_copy", "STAT1 (a)", "TP53 (a)", "B cell", "T cell CD4+", "T cell CD8+",
             "Macrophage", "NK cell", "Cancer associated fibroblast", "Endothelial cell",
             "uncharacterized cell", "STAT1", "MYC", "IFNG", "CIN70_Carter2006")
    pdf(args$plotfile, 8, 6)
    corrplot::corrplot(cor(x), tl.cex=0.5)
    gu$pcor(x) %>% gu$plot_pcor_net(node_size=4, edge_size=2.5)
    gu$plot_bootstrapped_pcor(x[,xsub], node_size=4)
    gu$plot_bootstrapped_pcor(x[,setdiff(xsub, "TP53 (a)")], node_size=4)
    gu$plot_bootstrapped_pcor(x[,sub("rev24", "rev48", xsub)], node_size=4)
    gu$plot_bootstrapped_pcor(x[,sub("rev24", "rev48", xsub) %>% setdiff("TP53 (a)")], node_size=4)
    print(p1)
    print(lm_plot(dset, aes(x=myc_copy, y=`Myc Targets V1`), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=`Interferon Gamma Response`), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=aneup_log2seg), covar="purity"))
    print(lm_plot(dset, aes(y=`Myc Targets V1`, x=aneup_abs), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=wt_ifn2_over_dmso), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=wt_rev24_over_dmso), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=wt_rev48_over_dmso), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=rev24_stat1_over_wt), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=rev48_stat1_over_wt), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=wt_ifn2_over_dmso), covar="rev24_stat1_over_wt"))
    print(lm_plot(dset, aes(x=aneup_abs, y=rev24_stat1_over_wt), covar="wt_ifn2_over_dmso")) # no cor
    print(lm_plot(dset, aes(x=wt_rev24_over_dmso, y=rev24_stat1_over_wt), covar="purity"))
    print(lm_plot(dset, aes(x=wt_rev48_over_dmso, y=rev24_stat1_over_wt), covar="Myc Targets V1"))
    print(lm_plot(dset, aes(x=wt_rev48_over_dmso, y=rev48_stat1_over_wt), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=rev48_stat1_over_wt), covar="CIN70_Carter2006"))
    print(lm_plot(dset, aes(x=CIN70_Carter2006, y=rev48_stat1_over_wt), covar="Myc Targets V1"))
    print(lm_plot(dset, aes(x=aneup_abs, y=wt_rev24_over_dmso), covar="purity")) # no cor
    print(lm_plot(dset, aes(x=aneup_log2seg, y=`Interferon Gamma Response`), covar="Myc Targets V1"))
    print(lm_plot(dset, aes(x=aneup_abs, y=`Interferon Gamma Response`), covar="Myc Targets V1"))
    print(lm_plot(dset, aes(x=rev24_stat1_over_wt, y=`purity`), covar="Interferon Gamma Response"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=purity), covar="Interferon Gamma Response"))
    print(lm_plot(dset, aes(x=purity, y=`Myc Targets V1`), covar="Interferon Gamma Response"))
    print(lm_plot(dset, aes(x=purity, y=`Myc Targets V1`), covar="rev24_stat1_over_wt"))
    dev.off()

    # split Myc targets low vs high (bimodal)
    # split IFN response low vs high (bimodal)

    library(survival)
    coxph(Surv(OS_time, vital_status) ~ purity+ aneup_log2seg,
          data=dset %>% filter(rev24_stat1_over_wt>0))# still sign after +IFNg resp/2h; <0 n.s. (and more surv)
    # still sign for p53_mut==0, not for !=0; both ns <0 (-> pt. that )
    coxph(Surv(OS_time, vital_status) ~ purity+ aneup_log2seg,
          data=dset %>% filter(rev24_stat1_over_wt>0, p53_mut==0))


    # aneup generally predictive of survival, MycV1 is not
    coxph(Surv(OS_time, vital_status) ~ age_at_diagnosis + aneup_log2seg, data=dset)
    coxph(Surv(OS_time, vital_status) ~ age_at_diagnosis + `Myc Targets V1`, data=dset)

    # can stratify with myc+stat1 act if p53 wt
    dset2 = dset %>%
        filter(p53_mut == 0) %>%
        mutate(class = case_when(
            `Myc Targets V1`< 0 ~ "nomyc", # <0.5 (high-myc) even has est<0 for stat1wt
            rev48_stat1_over_wt>0 ~ "stat1ko_myc",
            rev48_stat1_over_wt<0 ~ "stat1wt_myc"
        ))
    coxph(Surv(OS_time, vital_status) ~ age_at_diagnosis + class:aneup_log2seg, data=dset2) # sign w/ purity+, IFNg+, MycV1+; or all 3

    ggplot(dset2, aes(x=aneup_log2seg, y=purity)) + geom_point() + facet_wrap(~class) + geom_smooth(method="lm")
    ggplot(dset2, aes(x=class, y=purity, fill=cut(aneup_log2seg, c(-Inf,0.1,0.5,Inf)))) + geom_boxplot() #rm green
    ggplot(dset2, aes(x=class, y=purity, fill=cut(`Myc Targets V1`, c(-Inf,0,0.5,Inf)))) + geom_boxplot()
    ggplot(dset2 %>% filter(aneup_log2seg > 0.5), aes(x=class, y=purity)) + geom_boxplot()
    ggplot(dset2 %>% filter(CIN70_Carter2006 > 0.5), aes(x=class, y=purity)) + geom_boxplot()
    ggplot(dset2 %>% filter(`Myc Targets V1` > 0.5), aes(x=class, y=purity)) + geom_boxplot()
    dset2 %>% filter(`Myc Targets V1` > quantile(`Myc Targets V1`, 0.75)) %>% lm(purity ~ type+ class, data=.) %>% broom::tidy()

    # no diff in p53 mut
    dset3 = dset %>%
        filter(p53_mut == 1) %>%
        mutate(class = case_when(
            `Myc Targets V1`<0 ~ "nomyc", # <0.5 (high-myc) has stat1wt least surv by far
            rev48_stat1_over_wt>0 ~ "stat1ko_myc",
            rev48_stat1_over_wt<0 ~ "stat1wt_myc"
        ))
    coxph(Surv(OS_time, vital_status) ~ age_at_diagnosis + class:aneup_log2seg, data=dset3)



    dset4 = dset %>%
        filter(p53_mut == 0) %>%
        mutate(class = case_when(
            wt_rev48_over_dmso<0 & rev48_stat1_over_wt>0 ~ "CIN_stat1ko",
            wt_rev48_over_dmso>0 ~ "CIN",
            wt_rev48_over_dmso<0 & rev48_stat1_over_wt<0 ~ "noCIN"
        ), class=relevel(factor(class), "noCIN"))
    coxph(Surv(OS_time, vital_status) ~ age_at_diagnosis + class, data=dset4) # + :aneup_log2seg

    dset5 = dset %>%
        filter(p53_mut == 1) %>%
        mutate(class = case_when(
            wt_rev48_over_dmso<0 & rev48_stat1_over_wt>0 ~ "CIN_stat1ko",
            wt_rev48_over_dmso>0 ~ "CIN",
            wt_rev48_over_dmso<0 & rev48_stat1_over_wt<0 ~ "noCIN"
        ), class=relevel(factor(class), "noCIN"))
    coxph(Surv(OS_time, vital_status) ~ age_at_diagnosis + class, data=dset5) # + :aneup_log2seg

    bind_rows(list(p53wt=dset4, p53mut=dset5), .id="type") %>%
    mutate(MycV1 = cut(`Myc Targets V1`, c(-Inf,0,Inf))) %>%
    ggplot(aes(x=class, fill=class, y=1-purity)) +
        geom_boxplot(outlier.shape=NA) + facet_grid(MycV1 ~ type)

    ggplot(dset, aes(x=`Buccitelli_up`, y=rev48_stat1_over_wt)) + geom_point()
    lm_plot(dset, aes(x=aneup_log2seg, y=rev48_stat1_over_wt), covar="CIN70_Carter2006")

    #todo: non-purity STAT1 act? +"separate from Myc act"?
})
