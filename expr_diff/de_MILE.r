library(dplyr)
io = import('io')
sys = import('sys')
util = import('./util')
gset = import('data/genesets')

args = sys$cmd$parse(
    opt('e', 'eset', 'gene expression RData', 'eset_MILE.RData'),
    opt('f', 'config', 'yaml', '../config.yaml'),
    opt('o', 'outfile', 'results RData', 'de_MILE.RData'),
    opt('p', 'plotfile', 'pdf', 'de_MILE.pdf'),
    arg('sets', 'gene set .RData', arity='*',
        list.files("../data/genesets/mouse", "\\.RData", full.names=TRUE)))

dset = io$load(args$eset)
keep = !is.na(dset$meta$type)
expr = dset$expr[,keep]
#lineage = narray::mask(dset$meta$lineage, along=2) + 0
type = cbind('(Intercept)' = 1, narray::mask(dset$meta$type[keep]))
#aneuploidy = pmin(dset$meta$aneuploidy[keep], 0.25)
aneuploidy = (dset$meta$annot[keep] == "ALL with hyperdiploid karyotype") + 0

res = data.frame(gene_name = rownames(expr), mean_expr = rowMeans(expr)) %>%
    mutate(fit = purrr::map(gene_name, function(g)
        broom::tidy(lm(expr[g,] ~ type + aneuploidy)))) %>%
    tidyr::unnest() %>%
    filter(term != "(Intercept)") %>%
    group_by(term) %>%
    mutate(adj.p = 2^-abs(statistic)) %>%
    arrange(adj.p) %>%
    tidyr::nest()

res = setNames(res$data, res$term)

args$sets = sub("mouse/", "human/", args$sets)
sets = io$load(args$sets) %>%
    setNames(tools::file_path_sans_ext(basename(args$sets))) %>%
    lapply(function(x) gset$filter(x, min=5, valid=rownames(expr)))

hl = toupper(io$read_yaml(args$config)$highlight_de)

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
