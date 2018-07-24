b = import('base')
io = import('io')
sys = import('sys')
util = import('./qc_rnaseq')

args = sys$cmd$parse(
    opt('p', 'plotfile', 'pdf', 'assemble.pdf'),
    opt('o', 'outfile', 'RData', 'assemble.RData'))

load_file = function(fname) {
    cont = io$load(fname)
    batch = b$grep("batch([0-9]+)", fname)
    ee = cont$counts
    colnames(ee) = paste(batch, colnames(ee), sep="_")
    list(idx = cont$idx %>%
            mutate(sample = paste(batch, sample, sep="_"),
                   batch = batch),
         counts = ee)
}
cont = lapply(sprintf("Mad2+PB_batch%i.RData", 1:3), load_file)

idx = lapply(cont, function(x) x$idx) %>%
    dplyr::bind_rows()
counts = lapply(cont, function(x) x$counts) %>%
    narray::stack(along=2)
expr = util$rnaseq$vst(counts)
dd = as.matrix(dist(t(expr)))

annot = as.data.frame(idx)
rownames(annot) = idx$sample
annot$hist_nr = NULL
annot$sample = NULL

pdf(args$plotfile, 16, 14)
pheatmap::pheatmap(dd, col=colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
                   annotation=annot)
print(util$plot_pca(expr, idx))
print(util$plot_tsne(expr, idx))
dev.off()

save(idx, counts, expr, file=args$outfile)
