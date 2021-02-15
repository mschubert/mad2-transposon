library(dplyr)
sys = import('sys')
plt = import('plot')
ma = import('process/microarray')

plot_pca = function(meta, expr, color) {
    pca = prcomp(t(expr), scale=FALSE)
    ggplot(cbind(meta, pca$x), aes_string(x="PC1", y="PC2", color=color, shape="series")) +
        geom_point(size=5) +
        ggrepel::geom_text_repel(aes_string(label=color), color="black") +
        labs(x = sprintf("PC1 (%.1f%%)", summary(pca)$importance[2,1]*100),
             y = sprintf("PC2 (%.1f%%)", summary(pca)$importance[2,2]*100),
             title = "PCA plot (linear)")
}

args = sys$cmd$parse(
    opt('a', 'annot', 'yaml', 'dset.yaml'),
    opt('o', 'outfile', 'rds', 'dset.rds'),
    opt('p', 'plotfile', 'pdf', 'dset.pdf')
)

annot = yaml::read_yaml(args$annot)
adf = stack(annot$cells)

# Load mouse hematopoeitic expression with good resolution
# http://servers.binf.ku.dk/bloodspot/?gene=RUNX1&dataset=nl_mouse_data
exps = c("E-GEOD-14833", "E-GEOD-6506") %>%
    lapply(ArrayExpress::ArrayExpress) %>%
    ma$qc() %>%
#    ma$normalize() %>% # AE package already does this
    ma$annotate(summarize="ensembl_gene_id", dset="mmusculus_gene_ensembl")

expr = lapply(exps, as.matrix) %>%
    narray::stack(along=2)

meta = sapply(exps, Biobase::pData, simplify=FALSE) %>%
    bind_rows(.id = "series") %>%
    transmute(series = series,
              fname = Array.Data.File,
              cells = ifelse(is.na(Comment..Sample_source_name.),
                           Scan.Name,
                           Comment..Sample_source_name.)) %>%
    mutate(cells = sub("[ _][12]$", "", cells),
           cells = sub(" Mouse bone marrow", "", cells, fixed=TRUE),
           cells = sub("Long-term hematopoietic stem cells ", "", cells, fixed=TRUE),
           annot = adf$ind[match(cells, adf$values)])

narray::intersect(expr, meta$fname, along=2)
expr = sva::ComBat(expr, batch=meta$series, par.prior=TRUE)

p1 = plot_pca(meta, expr, "cells")
meta = meta %>% filter(!is.na(annot))
narray::intersect(expr, meta$fname, along=2)
p2 = plot_pca(meta, expr, "annot")

pdf(args$plotfile, 15, 10)
print(p1)
print(p2)
dev.off()

saveRDS(list(meta=meta, expr=expr), file=args$outfile)
