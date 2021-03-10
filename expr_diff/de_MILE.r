library(dplyr)
sys = import('sys')
util = import('./util')
gset = import('data/genesets')
viper = import('../expr_cor/viper')

args = sys$cmd$parse(
    opt('e', 'eset', 'gene expression rds', 'eset_MILE.rds'),
    opt('f', 'config', 'yaml', '../config.yaml'),
    opt('n', 'network', 'rds', '../data/networks/E-GEOD-13159.rds'),
    opt('o', 'outfile', 'results rds', 'de_MILE.rds'),
    opt('p', 'plotfile', 'pdf', 'de_MILE.pdf'),
    arg('sets', 'gene set .rds', arity='*',
        list.files("../data/genesets/mouse", "\\.rds", full.names=TRUE))
)

dset = readRDS(args$eset)
keep = !is.na(dset$meta$type)
expr = dset$expr[,keep]
type = cbind(narray::mask(dset$meta$type[keep]),
     Hyperdip = (dset$meta$annot[keep] == "ALL with hyperdiploid karyotype")) + 0
type[,"Hyperdip"][type[,"B_like"] == 0] = NA
aneuploidy = pmin(dset$meta$aneuploidy[keep], 0.25)
net = readRDS(args$network)
tf_net = filter(net, Target %in% Regulator)

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
sets = readRDS(args$sets) %>%
    setNames(tools::file_path_sans_ext(basename(args$sets))) %>%
    lapply(function(x) gset$filter(x, min=5, valid=rownames(expr)))

hl = toupper(yaml::read_yaml(args$config)$highlight_de)

pdf(args$plotfile)
for (rname in names(res)) {
    message(rname)
    print(util$plot_volcano(res[[rname]], hl) + ggtitle(rname))

    if (grepl("aneuploidy_in_", rname)) {
        tname = sub("aneuploidy_in_", "", rname)
        cmp = ifelse(type[,tname], aneuploidy, NA)
        aq = quantile(cmp, c(0, 0.33, 0.66, 1), na.rm=TRUE)
        ac = cut(cmp, aq, include.lowest=TRUE)
        levels(ac) = c("low", NA, "high")
        cmp = as.integer(ac == "high")
    } else
        cmp = type[,rname]
    dviper = viper$diff_viper(expr, net, cmp, fdr=0.3) %>% head(100)
    dcor = viper$diff_cor(expr, tf_net, cmp, fdr=0.3)
    print(viper$plot_subnet(dviper, dcor, highlight=hl) + ggtitle(rname))

    for (sname in names(sets)) {
        title = paste(rname, sname)
        message(title)
        print(util$plot_gset(res[[rname]], sets[[sname]]) + ggtitle(title))
    }
}
dev.off()

saveRDS(res, file=args$outfile)
