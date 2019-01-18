library(dplyr)
io = import('io')
sys = import('sys')
util = import('./util')
gset = import('data/genesets')

args = sys$cmd$parse(
    opt('e', 'eset', 'gene expression RData', 'eset_Mad2PB+multiCIS.RData'),
    opt('f', 'config', 'yaml', '../config.yaml'),
    opt('n', 'network', 'RData', '../data/networks/E-GEOD-13159.RData'), # ignored
    opt('o', 'outfile', 'results RData', 'de_Mad2PB+multiCIS.RData'),
    opt('p', 'plotfile', 'pdf', 'de_Mad2PB+multiCIS.pdf'),
    arg('sets', 'gene set .RData', arity='*',
        list.files("../data/genesets/mouse", "\\.RData", full.names=TRUE)))

dset = io$load(args$eset)
eset = dset$eset
fml = paste("~", paste(c("Tcell", "Other", dset$inc_ins), collapse=" + "))
res = util$do_wald(eset, as.formula(fml))

sets = io$load(args$sets) %>%
    setNames(tools::file_path_sans_ext(basename(args$sets))) %>%
    lapply(function(x) gset$filter(x, min=5, valid=rownames(eset)))

hl = io$read_yaml(args$config)$highlight_de

pdf(args$plotfile)
for (rname in names(res)) {
    message(rname)
    print(util$plot_volcano(res[[rname]], hl) + ggtitle(rname))
    for (sname in names(sets)) {
        title = paste(rname, sname)
        message(title)
        print(util$plot_gset(res[[rname]], sets[[sname]]) + ggtitle(title))
    }
}
dev.off()

save(res, file=args$outfile)
