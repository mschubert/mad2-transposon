library(dplyr)
sys = import('sys')
plt = import('plot')
util = import('./plot')

ccs = function(df) {
    df %>%
        filter(!is.na(cl), debris_gate) %>%
        util$cluster_centers(dset[[1]]$trans, inverse=FALSE)
}

args = sys$cmd$parse(
    opt('m', 'meta', 'tsv', '../meta/meta.tsv'),
    opt('o', 'outfile', 'rds', 'summarize.rds'),
    opt('p', 'plotfile', 'pdf', 'summarize.pdf'),
    arg('infiles', 'rds', arity='*', list.files(pattern="\\.rds"))
)

dset = lapply(args$infiles, readRDS)

res = dset[[1]]$res %>%
    lapply(ccs) %>%
    bind_rows(.id = "sample") %>%
    filter(grepl("^[0-9]{3}", sample))

annot = res %>%
    select(sample, cl, pct) %>%
    mutate(hist_nr = as.integer(sub("^([^ ]+).*", "\\1", sample)),
           sample = paste0(hist_nr, "s"),
           label = sprintf("%s.%s", sample, cl),
           clust = factor(clust$membership)) %>%
    left_join(readr::read_tsv(args$meta), by=c("sample", "hist_nr"))
#todo: fix type assignment with "t","s" samples

dmat = data.matrix(res[names(dset[[1]]$trans)]) #todo: do we want scatter here?
colnames(dmat) = dset[[1]]$meta$desc[match(colnames(dmat), dset[[1]]$meta$name)]
clust = igraph::cluster_louvain(scran::buildSNNGraph(t(dmat), k=4))
pca = prcomp(dmat, scale.=TRUE) # does this still center?
#umap2 = ... #todo

p1 = plt$pca(pca, aes(x=PC1, y=PC2), annot=annot, biplot=TRUE) +
    geom_point(aes(color=type, shape=clust, size=pct), alpha=0.8) +
    scale_size_area() +
    ggrepel::geom_text_repel(aes(label=label), size=2, max.overlaps=Inf) +
    theme_classic()

both = cbind(annot, dmat) %>%
    tidyr::gather("marker", "intensity", -(sample:type))
p2 = ggplot(both, aes(x=clust, y=intensity, group=clust)) +
    geom_point(aes(color=type, size=pct), alpha=0.3) +
    facet_wrap(~ marker, scales="free_y")

pdf(args$plotfile, 10, 9)
print(p1)
print(p2)
dev.off()
