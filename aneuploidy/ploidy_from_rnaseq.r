library(cowplot)
library(ggrepel)
io = import('io')
df = import('data_frame')
idmap = import('process/idmap')

# reads per chromosome for T-ALLs
dset = io$load('../data/rnaseq/assemble.RData')
chrs = idmap$gene(rownames(dset$counts), to="chromosome_name", dset="mmusculus_gene_ensembl")
reads = dset$counts[!is.na(chrs),]
chrs = chrs[!is.na(chrs)]
reads = narray::map(reads, along=1, sum, subsets=chrs)[as.character(1:19),]

# reads per chromosome for euploid controls
ctl = readr::read_tsv("../../aneuploidy/data/rnaseq/T-ALL_read_count.txt")
ctl_reads = data.matrix(ctl[,c("eT_p0","eT_p2")])
chrs = idmap$gene(ctl$ensembl_gene_id, to="chromosome_name", dset="mmusculus_gene_ensembl")
ctl_reads = ctl_reads[!is.na(chrs),]
chrs = chrs[!is.na(chrs)]
ctl_reads = narray::map(ctl_reads, along=1, sum, subsets=chrs)[as.character(1:19),]

# infer ploidy from read count difference
ploidy_multiples = reads / narray::crep(rowSums(ctl_reads), ncol(reads))
medians = narray::map(ploidy_multiples, along=1, median)
ploidy = 2 * ploidy_multiples / narray::rrep(medians, nrow(ploidy_multiples))
aneuploidy = colSums(abs(2 - ploidy)) #TODO: scale chrom length

# compare to measured ploidy using scWGS
scWGS2aneuploidy = function(sample_id) {
    mod = io$load(paste0("../../aneuploidy/data/singlecell_wgs/T-ALL/",
                         sample_id, ".RData"))
    score = AneuFinder::karyotypeMeasures(mod)$per.chromosome
    setNames(score[,"Aneuploidy"], rownames(score))
}
wgs = c("T401", "T419", "S413") %>%
    lapply(scWGS2aneuploidy) %>%
    setNames(c("401t", "419t", "413s")) %>%
    narray::stack(along=2)

names(dimnames(ploidy)) = names(dimnames(wgs)) = c("chr", "sample")
df_inf = reshape2::melt(ploidy) %>%
    mutate(rna = value) %>%
    select(-value)
df_obs = reshape2::melt(wgs) %>%
    mutate(wgs = value) %>%
    select(-value)
dfs = inner_join(df_inf, df_obs, by=c("chr", "sample")) %>%
    mutate(label = paste(sample, chr, sep=":"))
mod = lm(rna ~ wgs, data=dfs)
p = ggplot(dfs, aes(x=wgs, y=rna)) +
    geom_smooth(method="lm") +
    geom_point() +
    geom_text_repel(aes(label=label), size=3) +
    labs(x = "Ploidy from single-cell WGS data + AneuFinder",
         y = "Ploidy inferred from RNA-seq",
         title = sprintf("RNA-seq for ploidy inference (p=%.2g, r^2=%.2f)",
             mod %>% broom::tidy() %>% filter(term == "wgs") %>% pull(p.value),
             mod %>% broom::glance() %>% pull(r.squared)))

pdf("ploidy_from_rnaseq.pdf")
print(p)
dev.off()

save(ploidy, aneuploidy, wgs, file="ploidy_from_rnaseq.RData")
