library(dplyr)
sys = import('sys')
plt = import('plot')
ma = import('process/microarray')

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'dset.rds'),
    opt('p', 'plotfile', 'pdf', 'dset.pdf')
)

# Load mouse hematopoeitic expression with good resolution
# http://servers.binf.ku.dk/bloodspot/?gene=RUNX1&dataset=nl_mouse_data
exps = c("E-GEOD-14833", "E-GEOD-6506") %>%
    lapply(ArrayExpress::ArrayExpress) %>%
    ma$qc() %>%
#    ma$normalize() %>% # AE package already does this
    ma$annotate(summarize="ensembl_gene_id")

expr = lapply(exps, as.matrix) %>%
    narray::stack(along=2)

meta = sapply(exps, Biobase::pData, simplify=FALSE) %>%
    bind_rows(.id = "series") %>%
    transmute(series = series,
              fname = Array.Data.File,
              cells = ifelse(is.na(Comment..Sample_source_name.),
                           Scan.Name,
                           Comment..Sample_source_name.))

narray::intersect(expr, meta$fname, along=2)

expr = sva::ComBat(expr, batch=meta$series, par.prior=TRUE)

pca = prcomp(t(expr), scale=FALSE)
p1 = ggplot(cbind(meta, pca$x), aes(x=PC1, y=PC2, color=cells, shape=series)) +
    geom_point(size=5) +
    geom_text_repel(aes(label=cells), color="black") +
    labs(x = sprintf("PC1 (%.1f%%)", summary(pca)$importance[2,1]*100),
         y = sprintf("PC2 (%.1f%%)", summary(pca)$importance[2,2]*100),
         title = "PCA plot (linear)")

pdf(args$plotfile)
print(p1)
dev.off()

saveRDS(list(meta=meta, expr=expr), file=args$outfile)
