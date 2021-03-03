library(dplyr)
sys = import('sys')
plt = import('plot')
util = import('./plot')

ccs = function(df) {
    df %>%
        filter(!is.na(cl), debris_gate) %>%
        util$cluster_centers(trans, inverse=FALSE)
}

args = sys$cmd$parse(
    opt('m', 'meta', 'tsv', '../meta/meta.tsv'),
    opt('o', 'outfile', 'rds', 'summarize.rds'),
    opt('p', 'plotfile', 'pdf', 'summarize.pdf'),
    arg('infiles', 'rds', arity='*', list.files(pattern="\\.rds"))
)

dset = lapply(args$infiles, readRDS)

meta = dset[[1]]$meta # same
trans = dset[[1]]$trans # same
res = c(dset[[1]]$res, dset[[2]]$res) %>%
    lapply(ccs) %>%
    bind_rows(.id = "sample") %>%
    filter(grepl("^[0-9]{3}", sample))

dmat = data.matrix(res[names(trans)]) #todo: do we want scatter here?
colnames(dmat) = meta$desc[match(colnames(dmat), meta$name)]
clust = igraph::cluster_louvain(scran::buildSNNGraph(t(dmat), k=10))
pca = prcomp(dmat, scale.=TRUE) # does this still center?
umap2 = uwot::umap(dmat, n_components=2)
colnames(umap2) = c("umap1", "umap2")

annot = res %>%
    select(sample, cl, pct) %>%
    mutate(hist_nr = as.integer(sub("^([^ ]+).*", "\\1", sample)),
           sample = paste0(hist_nr, "s"),
           label = sprintf("%s.%s", sample, cl),
           clust = factor(clust$membership)) %>%
    left_join(readr::read_tsv(args$meta), by=c("sample", "hist_nr"))
#todo: fix type assignment with "t","s" samples

p1 = plt$pca(pca, aes(x=PC1, y=PC2), annot=annot, biplot=TRUE) +
    geom_point(aes(color=type, shape=clust, size=pct), alpha=0.8) +
    scale_size_area() +
    ggrepel::geom_text_repel(aes(label=label), size=2, max.overlaps=Inf) +
    theme_classic()

p2 = ggplot(cbind(annot, umap2), aes(x=umap1, y=umap2)) +
    geom_point(aes(color=type, shape=clust, size=pct), alpha=0.8) +
    scale_size_area() +
    ggrepel::geom_text_repel(aes(label=label), size=2, max.overlaps=Inf) +
    theme_classic()

both = cbind(annot, dmat) %>%
    tidyr::gather("marker", "intensity", -(sample:type))
p3 = ggplot(both, aes(x=clust, y=intensity, group=clust)) +
    geom_point(aes(color=type, size=pct), alpha=0.3) +
    facet_wrap(~ marker, scales="free_y")

p4 = ggplot(both, aes(x=intensity)) +
    geom_density() +
    facet_wrap(~ marker)

pdf(args$plotfile, 10, 9)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()
