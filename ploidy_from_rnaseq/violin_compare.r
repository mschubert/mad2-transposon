library(reshape2)
library(dplyr)
library(ggplot2)
io = import('io')
idmap = import('process/idmap')

load('../data/rnaseq/assemble.RData') # counts, expr, idx
ctl = edgeR::cpm(io$read_table("../data/rnaseq/read_count.txt"))[,c("eT_p0", "eT_p2")]
samples = edgeR::cpm(counts)
narray::intersect(samples, ctl, along=1)
keep = rowMeans(ctl) > 1 & rowMeans(ctl) < 100 &
       rowMeans(samples) > 1 & rowMeans(samples) < 100
scaled = 2*(samples / rowMeans(ctl))[keep,]

df = melt(scaled) %>%
	mutate(Var1 = as.character(Var1),
           Var2 = as.character(Var2))
df$chr = idmap$gene(df$Var1, from="ensembl_gene_id", to="chromosome_name", dset="mmusculus_gene_ensembl")
df = df %>%
    filter(chr %in% c(1:19,'X')) %>% # no Y in ctl
    mutate(chr = factor(chr, levels=c(1:19,'X')))

p = ggplot(df, aes(x=chr, y=value)) +
    geom_violin() +
    geom_hline(yintercept=c(1:8), linetype="dashed", alpha=0.3) +
    stat_summary(fun.y="median", color="red", geom="point") +
    facet_grid(Var2~chr, scales="free_x") +
    scale_y_continuous(name="ploidy", trans="log", breaks=1:8) +
    coord_cartesian(ylim=c(1, 4)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

pdf("violin_compare.pdf", 10,40)
print(p)
dev.off()
