library(dplyr)
sys = import('sys')
util = import('../../expr_diff/util')

args = sys$cmd$parse(
    opt('d', 'dset', 'rds', 'dset.rds'),
    opt('o', 'outfile', 'rds', 'diff_expr.rds'),
    opt('p', 'plotfile', 'pdf', 'diff_expr.pdf'))

dset = readRDS(args$dset)
eset = dset$eset %>%
    DESeq2::estimateSizeFactors()
samples = dset$index

genotype = eset[,eset$time == 0]
res = util$do_wald(genotype, ~ genotype, ex="genotype")

over_dmso = function(genotype, time, treatment) {
    cur = eset[,eset$time %in% time & eset$genotype == genotype & eset$treatment %in% c(treatment, "dmso")]
    colData(cur) = droplevels(colData(cur))
    res = util$do_wald(cur, ~ treatment, ex="treatment")
}
res$wt_ifn2_over_dmso = over_dmso("wt", c("0", "2"), "ifng")
res$wt_rev24_over_dmso = over_dmso("wt", "24", "rev")
res$wt_rev48_over_dmso = over_dmso("wt", "48", "rev")

res$stat1_rev24_over_dmso = over_dmso("stat1", "24", "rev")
res$stat1_rev48_over_dmso = over_dmso("stat1", "48", "rev")

# rep 1+2 too different, only DE genes if replicate is regressed out
cgas_over_dmso = function(genotype, time, treatment) {
    cur = eset[,eset$time %in% time & eset$genotype == genotype & eset$treatment %in% c(treatment, "dmso")]
    colData(cur) = droplevels(colData(cur))
    res = util$do_wald(cur, ~ replicate + treatment, ex="treatment")
}
res$cgas_rev24_over_dmso = cgas_over_dmso("cgas", "24", "rev")

over_wt = function(genotype, time, treatment) {
    cur = eset[,eset$time %in% time & eset$genotype %in% genotype & eset$treatment %in% treatment]
    colData(cur) = droplevels(colData(cur))
    res = util$do_wald(cur, ~ genotype, ex="genotype")
}
res$rev24_stat1_over_wt = over_wt(c("wt", "stat1"), "24", "rev")
res$rev48_stat1_over_wt = over_wt(c("wt", "stat1"), "48", "rev")

#TODO: gene sets

hl = c()

pdf(args$plotfile)
for (rname in names(res)) {
    message(rname)
    print(util$plot_volcano(res[[rname]], hl) + ggtitle(rname))
}
dev.off()
