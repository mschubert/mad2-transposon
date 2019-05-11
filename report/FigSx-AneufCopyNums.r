library(dplyr)
library(cowplot)
io = import('io')
seq = import('seq')
anp = import('seq/aneuploidy')

aset = io$load("../ploidy_compare/analysis_set.RData")
mrg = io$read_table("../ploidy_compare/analysis_set_merge.tsv", header=TRUE)
m1 = io$load("../data/wgs/30cellseq_batch1.RData")
m2 = io$load("../data/wgs/30cellseq_batch2.RData")
models = utils::modifyList(m1, m2)
asm = GenomeInfoDb::seqinfo(m1[[1]]$segments)

adj_segs = io$load("../data/wgs/30cellseq.RData")$segments
gr = split(adj_segs, adj_segs$sample) 

#FIXME: this needs go work for seq$aneuploidy
#aneups = anp$aneuploidy(adj_segs, assembly=asm, sample="sample")
aneups = adj_segs %>%
    filter(seqnames != "X") %>%
    group_by(sample) %>%
    summarize(aneup = weighted.mean(abs(2-ploidy), w=width)) %>%
    arrange(aneup) %>%
    filter(sample %in% c(aset$sample, mrg$subset)) %>%
    mutate(sample = factor(sample, levels=sample))

models = models[levels(aneups$sample)]
for (i in seq_along(models)) {
    int = as.integer(round(gr[[models[[i]]$ID]]$ploidy))
    cn = paste0(int, "-somy")
    cn = factor(cn, levels=levels(models[[i]]$segments$state))
    models[[i]]$segments$state = cn
    models[[i]]$segments$copy.number = int
}

pdf("FigSx-AneufCopyNums.pdf", 18, 12)
AneuFinder::heatmapGenomewide(models, cluster=FALSE)
dev.off()
