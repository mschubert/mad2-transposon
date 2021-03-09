library(dplyr)
library(ggplot2)
sys = import('sys')

args = sys$cmd$parse(
    opt('d', 'dset', 'rds', 'dset.rds'),
    opt('g', 'go', 'rds', '../data/gsva/tcga-brca/GO_Biological_Process_2020.rds'),
    opt('p', 'plotfile', 'pdf', 'explain_myc_pur_cor.pdf')
)

ds = readRDS(args$dset)
dset = cbind(ds$meta[c("purity")], as.data.frame(ds$dmat))
go = t(readRDS(args$go))
narray::intersect(dset, go, along=1)

both = as.data.frame(cbind(dset, go))


test_one = function(covar) {
    full = broom::tidy(lm(purity ~ covar + `Myc Targets V2`, data=both)) %>%
        filter(term != "(Intercept)") %>%
        mutate(term = sub("`Myc Targets V2`", "Myc", term)) %>%
        select(term, est=estimate, stat=statistic, p.value) %>%
        tidyr::pivot_wider(names_from=c(term), values_from=c(estimate:p.value)) %>%
        mutate(delta_stat = red_stat - stat_Myc)
}
naive = broom::tidy(lm(purity ~ `Myc Targets V2`, data=both))
red_stat = naive %>% filter(term == "`Myc Targets V2`") %>% pull(statistic)

res = clustermq::Q(test_one, covar=narray::split(both[,-1], along=2), n_jobs=0,
                   export = list(naive=naive, both=both, red_stat=red_stat),
                   pkgs = c("dplyr")) %>%
    setNames(colnames(both)[-1]) %>%
    bind_rows(.id="covar") %>%
    arrange(delta_stat)

res %>% filter(estimate_covar > 0) %>% select(covar, delta_stat)
res %>% filter(estimate_covar > 0, !grepl("^GO:", covar))
