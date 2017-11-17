library(dplyr)
library(factoextra)
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
ggplot(cbind(samples, pca$x), aes(x=PC1, y=PC2, color=tissue, shape=`Early/ Late`)) +
    geom_point(size=5) +
    geom_text(aes(label=id), color="black") +
    labs(x = sprintf("PC1 (%.1f%%)", summary(pca)$importance[2,1]*100),
         y = sprintf("PC2 (%.1f%%)", summary(pca)$importance[2,2]*100))

tsne = Rtsne::Rtsne(t(expr), perplexity=10)
cbind(samples, x=tsne$Y[,1], y=tsne$Y[,2]) %>%
    ggplot(aes(x=x, y=y, color=tissue, shape=`Early/ Late`)) +
    geom_point(size=5) +
    geom_text(aes(label=id), color="black")

# save in merged object
