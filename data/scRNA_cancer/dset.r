library(dplyr)
library(ggplot2)
gset = import('genesets')

plot_gset = function(scde, field) {
    sf = rlang::sym(field)
    res = scde %>%
        select(dset, !! sf) %>%
        tidyr::unnest(!! sf)

    sums = res %>% group_by(label) %>%
        summarize(mod = list(broom::tidy(lm(statistic ~ 1)))) %>%
        tidyr::unnest(mod) %>%
        arrange(estimate)

    res$label = factor(res$label, levels=sums$label)

    mat = narray::construct(estimate ~ dset + label, data=res)
    mat[is.na(mat)] = 0
    clust = uwot::umap(mat, n_components=1, n_neighbors=5)
    res$dset = factor(res$dset, levels=c(rownames(mat)[order(clust)]))

    ggplot(res, aes(x=label, y=dset)) +
        geom_raster(aes(fill=statistic)) +
        scale_fill_fermenter(palette="RdBu", breaks=c(-Inf,-5,-3,-1.5,-0.5,0.5,1.5,3,5,Inf)) +
        theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5))
}

bpath = "/home/p.batistatan/main/NKI-CIN-internship-main-backup"

dsets = c(
    "BakhoumEtAl2018",
#    "ChungEtAl2017", # unclear CIN split
    "DarmanisEtAl2017",
    "GiustacchiniEtAl2017",
    "KaraayvazEtAl2018",
    "LeeEtAl2020",
    "NelsonEtAl2020",
    "PuramEtAl2017",
    "SunEtAl2020",
#    "TijhuisEtAl202X",
    "TiroshEtAl2016"
)

diff_expr = file.path(bpath, dsets, "output_DE/df.deseq.results.pseudobulkWald.csv") %>%
    lapply(read.table, sep=";", dec=",", header=TRUE) %>%
    lapply(as_tibble) %>%
    lapply(. %>% filter(baseMean >= 100))

sets = list(
    MSigDB_Hallmark_2020 = gset$get_human("MSigDB_Hallmark_2020"),
    DoRothEA = gset$get_human("DoRothEA", conf="A"),
    BT549 = readRDS("../genesets/human/stat1_ko.rds")
)

res = tibble(dset = dsets, genes = diff_expr) %>%
    rowwise() %>%
    mutate(MSigDB_Hallmark_2020 = list(gset$test_lm(genes, sets[[1]])),
           DoRothEA = list(gset$test_lm(genes, sets[[2]])),
           BT549 = list(gset$test_lm(genes, sets[[3]])))

saveRDS(res, file="dset.rds")

plots = lapply(names(sets), function(s) plot_gset(res, s))
pdf("dset.pdf", 12, 6)
for (p in plots)
    print(p)
dev.off()
