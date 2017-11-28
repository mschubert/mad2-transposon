library(cowplot)
library(ggridges)
io = import('io')
idmap = import('process/idmap')

# load gene expression
dset = io$load('../data/rnaseq/assemble.RData')
expr = dset$expr[rowSums(dset$counts) >= 10,]

# map it to chromosomes
rownames(expr) = idmap$gene(rownames(expr), to="chromosome_name", dset="mmusculus_gene_ensembl")
expr = expr[rownames(expr) %in% 1:19,]
names(dimnames(expr)) = c("chromosome", "sample")
chrs = reshape2::melt(expr) %>%
    split(.$chromosome)

# load known ploidies
# ...

# per chromosome: ggridges across samples
plot_fun = function(df) {
    # order by median value
    ord = df %>%
        group_by(sample) %>%
        summarize(median = median(value)) %>%
        arrange(median) #TOOD: merge known ploidies, z-transform for 1=ploidy

    # plot from low to high
    df %>%
        mutate(sample = factor(sample, levels=rev(ord$sample))) %>%
        ggplot(aes(x=value, y=sample)) +
            geom_density_ridges() +
            geom_point(data=ord, aes(x=median, y=sample), shape=5)

    #TOOD: plot inferred ploidies, align using cowplot
}

plots = lapply(chrs, plot_fun)

pdf("ridges.pdf", width=8, height=6)
for (i in seq_along(plots))
    print(plots[[i]] + ggtitle(paste("chromosome", names(plots)[i])))
dev.off()
