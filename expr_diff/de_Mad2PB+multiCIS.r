library(dplyr)
sys = import('sys')
util = import('./util')
gset = import('data/genesets')

args = sys$cmd$parse(
    opt('e', 'eset', 'gene expression rds', 'eset_Mad2PB+multiCIS.rds'),
    opt('f', 'config', 'yaml', '../config.yaml'),
    opt('n', 'network', 'rds', '../data/networks/E-GEOD-13159.rds'), # ignored
    opt('o', 'outfile', 'results rds', 'de_Mad2PB+multiCIS.rds'),
    opt('p', 'plotfile', 'pdf', 'de_Mad2PB+multiCIS.pdf'),
    arg('sets', 'gene set .rds', arity='*',
        list.files("../data/genesets/mouse", "\\.rds", full.names=TRUE))
)

dset = readRDS(args$eset)
eset = dset$eset
fml = paste("~", paste(c("Tcell", "Other", dset$inc_ins), collapse=" + "))
res = util$do_wald(eset, as.formula(fml))

sets = readRDS(args$sets) %>%
    setNames(tools::file_path_sans_ext(basename(args$sets))) %>%
    lapply(function(x) gset$filter(x, min=5, valid=rownames(eset)))

hl = yaml::read_yaml(args$config)$highlight_de

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

saveRDS(res, file=args$outfile)
