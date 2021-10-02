# data frame with:
# * myc copy number ,expression level, target expressioin
# * immune cells
# * subtype
# * survival
library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')

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
    ds$meta$OS_time = ifelse(ds$meta$OS_time > 365*5, NA, ds$meta$OS_time)
    dset = cbind(ds$meta, as.data.frame(ds$dmat))

    dens = dset %>%
        select(Sample, purity, aneuploidy, `Aneuploidy Score`, `Interferon Gamma Response`,
               `Myc Targets V1`, `Leukocyte Fraction`, `Stromal Fraction`, `Lymphocytes`,
               wt_rev24_over_dmso, wt_rev48_over_dmso, rev24_stat1_over_wt, rev48_stat1_over_wt) %>%
        tidyr::gather("field", "value", -Sample)

    p1 = ggplot(dens, aes(x=value)) +
        geom_density() +
        facet_wrap(~ field, scales="free")

    x = na.omit(data.matrix(dset[-1]))
    pdf(args$plotfile, 8, 6)
    corrplot::corrplot(cor(x))
    print(p1)
    print(lm_plot(dset, aes(x=myc_copy, y=`Myc Targets V1`), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=`Interferon Gamma Response`), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=`Aneuploidy Score`), covar="purity"))
    print(lm_plot(dset, aes(y=`Myc Targets V1`, x=aneuploidy), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=wt_ifn2_over_dmso), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=wt_rev24_over_dmso), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=wt_rev48_over_dmso), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=rev24_stat1_over_wt), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=rev48_stat1_over_wt), covar="purity"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=wt_ifn2_over_dmso), covar="rev24_stat1_over_wt"))
    print(lm_plot(dset, aes(x=aneuploidy, y=rev24_stat1_over_wt), covar="wt_ifn2_over_dmso")) # no cor
    print(lm_plot(dset, aes(x=wt_rev24_over_dmso, y=rev24_stat1_over_wt), covar="purity"))
    print(lm_plot(dset, aes(x=wt_rev48_over_dmso, y=rev24_stat1_over_wt), covar="Myc Targets V1"))
    print(lm_plot(dset, aes(x=wt_rev48_over_dmso, y=rev48_stat1_over_wt), covar="purity"))
    print(lm_plot(dset, aes(x=aneuploidy, y=wt_rev24_over_dmso), covar="purity")) # no cor
    print(lm_plot(dset, aes(x=`Aneuploidy Score`, y=`Interferon Gamma Response`), covar="Myc Targets V1"))
    print(lm_plot(dset, aes(x=aneuploidy, y=`Interferon Gamma Response`), covar="Myc Targets V1"))
    print(lm_plot(dset, aes(x=rev24_stat1_over_wt, y=`purity`), covar="Interferon Gamma Response"))
    print(lm_plot(dset, aes(x=`Myc Targets V1`, y=purity), covar="Interferon Gamma Response"))
    print(lm_plot(dset, aes(x=purity, y=`Myc Targets V1`), covar="Interferon Gamma Response"))
    print(lm_plot(dset, aes(x=purity, y=`Myc Targets V1`), covar="rev24_stat1_over_wt"))
    dev.off()

    # split Myc targets low vs high (bimodal)
    # split IFN response low vs high (bimodal)
})
