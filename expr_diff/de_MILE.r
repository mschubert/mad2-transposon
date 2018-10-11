library(dplyr)
io = import('io')
sys = import('sys')
util = import('./util')
gset = import('data/genesets')
viper = import('../link_cis_expr/cor_viper')

args = sys$cmd$parse(
    opt('e', 'eset', 'gene expression RData', 'eset_MILE.RData'),
    opt('f', 'config', 'yaml', '../config.yaml'),
    opt('n', 'network', 'RData', '../data/networks/E-GEOD-13159.RData'),
    opt('o', 'outfile', 'results RData', 'de_MILE.RData'),
    opt('p', 'plotfile', 'pdf', 'de_MILE.pdf'),
    arg('sets', 'gene set .RData', arity='*',
        list.files("../data/genesets/mouse", "\\.RData", full.names=TRUE)))

dset = io$load(args$eset)
keep = !is.na(dset$meta$type)
expr = dset$expr[,keep]
type = cbind(narray::mask(dset$meta$type[keep]),
     Hyperdip = (dset$meta$annot[keep] == "ALL with hyperdiploid karyotype")) + 0
type[,"Hyperdip"][type[,"B_like"] == 0] = NA
aneuploidy = pmin(dset$meta$aneuploidy[keep], 0.25)
aneup_tert = dset$meta[keep,] %>% # move this to a column in meta (@eset)?
    group_by(type) %>%
    mutate(tert = cut(aneuploidy, breaks=3, labels=FALSE)) %>%
    pull(tert) %>%
    factor()
levels(aneup_tert) = c("low", NA, "high")

edf = data.frame(gene_name = rownames(expr), mean_expr = rowMeans(expr), stringsAsFactors=FALSE)
ttf = expand.grid(gene_name = rownames(expr), term = colnames(type), stringsAsFactors=FALSE)
both = dplyr::inner_join(edf, ttf)
assocs_type = both %>%
    mutate(fit = purrr::map2(gene_name, term, function(g, t)
        r1 = broom::tidy(lm(expr[g,] ~ type[,t]))))
assocs_aneup = both %>%
    filter(term != "Hyperdip") %>%
    mutate(fit = purrr::map2(gene_name, term, function(g, t) {
            smps = !is.na(type[,t]) & type[,t] == 1
            broom::tidy(lm(expr[g,smps] ~ aneuploidy[smps]))
        }), term = paste0("aneuploidy_in_", term))
res = dplyr::bind_rows(assocs_type, assocs_aneup) %>%
    tidyr::unnest() %>%
    filter(term1 != "(Intercept)") %>%
    select(-term1) %>%
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

    if (grepl("aneuploidy_in_", rname))
        cmp = aneup_tert[type[,sub("aneuploidy_in_", "", rname)]]
    else
        cmp = type[,rname]
    dviper = viper$diff_viper(expr, net, tmat[,cond])
    dcor = viper$diff_cor(expr, tf_net, tmat[,cond])
    print(viper$plot_subnet(dviper, dcor, highlight=hl) +
          ggtitle(paste(rname, "aneup_high_vs_low_tertile")))

    for (sname in names(sets)) {
        title = paste(rname, sname)
        message(title)
        print(util$plot_gset(res[[rname]], sets[[sname]]) + ggtitle(title))
    }
}
dev.off()

save(res, file=args$outfile)
