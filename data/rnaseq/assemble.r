library(dplyr)
library(ggrepel)
plt = import('plot')
rnaseq = import('process/rna-seq')

# load samples and controls
samples = readr::read_tsv("samples.tsv") %>%
    mutate(Sample = id) %>%
    full_join(readr::read_tsv("groups.tsv"), by="Sample")
edf = readr::read_tsv("PB170929_raw_readcount_dedup_UMIs.txt")
expr = data.matrix(edf[,-1])
rownames(expr) = edf$Gene

expr = rnaseq$vst(expr)
narray::intersect(expr, samples$file, along=2)
colnames(expr) = samples$id

# do PCA/dim reduction plots to see how they cluster
pca = prcomp(t(expr), scale=FALSE)
p1 = ggplot(cbind(samples, pca$x), aes(x=PC1, y=PC2, color=`Tumour type`, shape=`Early/ Late`)) +
    geom_point(size=5) +
    geom_text_repel(aes(label=id), color="black") +
    labs(x = sprintf("PC1 (%.1f%%)", summary(pca)$importance[2,1]*100),
         y = sprintf("PC2 (%.1f%%)", summary(pca)$importance[2,2]*100),
         title = "PCA plot (linear)")

tsne = Rtsne::Rtsne(t(expr), perplexity=10)
p2 = cbind(samples, x=tsne$Y[,1], y=tsne$Y[,2]) %>%
    ggplot(aes(x=x, y=y, color=`Tumour type`, shape=`Early/ Late`)) +
    geom_point(size=5) +
    geom_text_repel(aes(label=id), color="black") +
    labs(x = "tsne 1",
         y = "tsne 2",
         title = "T-SNE plot (non-linear)")

pdf("overview.pdf", width=10, height=8)
print(p1)
print(p2)
dev.off()

# save in merged object
