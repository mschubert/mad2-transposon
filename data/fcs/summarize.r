library(dplyr)
sys = import('sys')
plt = import('plot')
util = import('./plot')

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'summarize.rds'),
    opt('p', 'plotfile', 'pdf', 'summarize.pdf'),
    arg('infiles', 'rds', arity='*', list.files(pattern="\\.rds"))
)

dset = lapply(args$infiles, readRDS)

#todo: transform coords
res = dset[[1]]$res %>%
    lapply(. %>% filter(!is.na(cl), debris_gate) %>% (util$cluster_centers)) %>% # todo: add fracs, pathology annots
    bind_rows(.id = "sample") %>%
    select(-Time, -debris_gate) %>%
    filter(grepl("^[0-9]{3}", sample))

dmat = data.matrix(res[-c(1:2)])
rownames(dmat) = paste(sub(" [0-9]{2}_[0-9]{3}", "", res$sample), res$cl, sep=".")

pca = prcomp(dmat, scale.=TRUE)
plt$pca(pca, aes(x=PC1, y=PC2), annot=data.frame(id=rownames(dmat)), biplot=TRUE) +
    geom_point() +
    geom_text(aes(label=id)) +
    theme_classic()
