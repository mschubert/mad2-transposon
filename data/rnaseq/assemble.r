b = import('base')
io = import('io')
util = import('./qc_rnaseq')

OUTFILE = "assemble.pdf"

load_file = function(fname) {
    cont = io$load(fname)
    batch = b$grep("batch([0-9]+)", fname)
    ee = cont$counts
    colnames(ee) = paste(batch, colnames(ee), sep="_")
    list(idx = idx %>% mutate(sample = paste(batch, sample, sep="_")),
         counts = ee)
}
cont = lapply(sprintf("Mad2+PB_batch%i.RData", 1:3), load_file)

idx = lapply(cont, function(x) x$idx) %>%
    dplyr::bind_rows()
expr = lapply(cont, function(x) x$counts) %>%
    narray::stack(along=2) %>%
    util$rnaseq$vst()
dd = as.matrix(dist(t(expr)))

pdf(OUTFILE, 14, 14)
pheatmap::pheatmap(dd, col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
print(util$plot_pca(expr, idx))
print(util$plot_tsne(expr, idx))
dev.off()
