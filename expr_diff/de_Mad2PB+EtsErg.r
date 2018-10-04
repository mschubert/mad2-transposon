library(dplyr)
io = import('io')
sys = import('sys')
util = import('./util')
gset = import('data/genesets')

args = sys$cmd$parse(
    opt('e', 'eset', 'gene expression RData', 'eset_Mad2PB+EtsErg.RData'),
    opt('c', 'cis', 'cis site RData', '../cis_analysis/poisson.RData'), # ignored (for now)
    opt('o', 'outfile', 'results RData', 'de_Mad2PB+EtsErg.RData'),
    opt('p', 'plotfile', 'pdf', 'de_Mad2PB+EtsErg.pdf'),
    arg('sets', 'gene set .RData', arity='*',
        list.files("../data/genesets", "\\.RData", full.names=TRUE)))

eset = io$load(args$eset)$eset

res = util$do_wald(eset, ~ group)
res$aneup0.3 = util$do_lrt(eset, ~ group + aneup0.3, ~ group)
res$`T+Ets_aneup` = util$do_wald(eset[,eset$group == "Ets1:thymus"], ~ aneup0.3)
res$`B+Ets_aneup` = util$do_wald(eset[,eset$group == "Ets1:spleen"], ~ aneup0.3)
res$`B+Erg_aneup` = util$do_wald(eset[,eset$group == "Erg:spleen"], ~ aneup0.3)

sets = io$load(args$sets) %>%
    setNames(tools::file_path_sans_ext(basename(args$sets))) %>%
    lapply(function(x) gset$filter(x, min=5, valid=rownames(eset)))

pdf(args$plotfile)
for (rname in names(res)) {
    message(rname)
    print(util$plot_volcano(res[[rname]]) + ggtitle(rname))
    for (sname in names(sets)) {
        title = paste(rname, sname)
        message(title)
        print(util$plot_gset(res[[rname]], sets[[sname]]) + ggtitle(title))
    }
}
dev.off()

save(res, file=args$outfile)
