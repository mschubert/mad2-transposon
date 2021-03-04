b = import('base')
sys = import('sys')
util = import('./qc_rnaseq')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('p', 'plotfile', 'pdf', 'assemble.pdf'),
    opt('o', 'outfile', 'rds', 'assemble.rds'))

load_file = function(fname) {
    cont = readRDS(fname)
    batch = b$grep("batch([0-9]+)", fname)
    ee = cont$counts
    colnames(ee) = paste(batch, colnames(ee), sep="_")
    list(idx = cont$idx %>%
            mutate(sample = paste(batch, sample, sep="_"),
                   batch = batch),
         counts = ee)
}
cont = lapply(sprintf("Mad2+PB_batch%i.rds", 1:3), load_file)

idx = lapply(cont, function(x) x$idx) %>%
    dplyr::bind_rows()
counts = lapply(cont, function(x) x$counts) %>%
    narray::stack(along=2)
expr = util$rnaseq$vst(counts)
dd = as.matrix(dist(t(expr)))

annot = as.data.frame(idx) %>% select(type, tissue, batch)
rownames(annot) = idx$sample

pdf(args$plotfile, 16, 14)
pheatmap::pheatmap(dd, col=colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255),
                   annotation=annot)
print(util$plot_pca(expr, idx))
print(util$plot_tsne(expr, idx))
dev.off()

colnames(counts) = sub("[0-9]_", "", colnames(counts))
colnames(expr) = sub("[0-9]_", "", colnames(expr))
keep = !duplicated(colnames(counts))
counts = counts[,keep]
expr = expr[,keep]
idx = idx[keep,]

genes = idmap$gene(rownames(expr), to="external_gene_name", dset="mmusculus_gene_ensembl")

saveRDS(list(idx=idx, counts=counts, expr=expr, genes=genes), file=args$outfile)
