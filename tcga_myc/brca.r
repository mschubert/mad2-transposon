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
        labs(subtitle = sprintf("p=%.2g (%.2g with %s)", m1$p.value, m2$p.value, covar))
}

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'dset.rds'),
    opt('o', 'outfile', 'rds', 'brca.rds'),
    opt('p', 'plotfile', 'pdf', 'brca.pdf')
)

ds = readRDS(args$infile)
dset = cbind(ds$meta, as.data.frame(ds$dmat))

lm_plot(dset, aes(x=`Myc Targets V2`, y=`Interferon Gamma Response`), covar="purity")
lm_plot(dset, aes(x=`Myc Targets V2`, y=aneup), covar="purity")
lm_plot(dset, aes(x=`Myc Targets V2`, y=purity), covar="Interferon Gamma Response")
dev.off()
