library(cowplot)
library(ggrepel)
io = import('io')
df = import('data_frame')
idmap = import('process/idmap')
aneufinder = import('tools/aneufinder')

# reads per chromosome for T-ALLs
dset = io$load('../../data/rnaseq/assemble.RData')
chrs = idmap$gene(rownames(dset$counts), to="chromosome_name", dset="mmusculus_gene_ensembl")
reads = dset$counts[!is.na(chrs),]
chrs = chrs[!is.na(chrs)]
reads = narray::map(reads, along=1, sum, subsets=chrs)[as.character(1:19),]
reads = reads / narray::rrep(colSums(reads), nrow(reads))

# reads per chromosome for euploid controls
ctl = readr::read_tsv("../../data/rnaseq/T-ALL_read_count.txt")
ctl_reads = data.matrix(ctl[,c("eT_p0","eT_p2")])
chrs = idmap$gene(ctl$ensembl_gene_id, to="chromosome_name", dset="mmusculus_gene_ensembl")
ctl_reads = ctl_reads[!is.na(chrs),]
chrs = chrs[!is.na(chrs)]
ctl_reads = narray::map(ctl_reads, along=1, sum, subsets=chrs)[as.character(1:19),]
ctl_reads = ctl_reads / narray::rrep(colSums(ctl_reads), nrow(ctl_reads))

# # compute ploidy z-scores
# # this does not work well because only 2 control samples (and non-cancer)
# z_mean = narray::map(ctl_reads, along=2, mean) %>% narray::crep(ncol(reads))
# z_sd = narray::map(ctl_reads, along=2, sd) %>% narray::crep(ncol(reads))
# z_ploidy = (reads - z_mean) / z_sd

# infer ploidy from read count difference
ploidy = 2 * reads / narray::crep(rowMeans(ctl_reads), ncol(reads))
aneuploidy = colSums(abs(2 - ploidy)) #TODO: scale chrom length

# compare to measured ploidy using scWGS
wgs = c("T401", "T419", "S413") %>%
    paste0("../../../dosage/data/singlecell_wgs/T-ALL/", ., ".RData") %>%
    io$load() %>%
    lapply(aneufinder$consensus_ploidy) %>%
    setNames(c("401t", "419t", "413s")) %>%
    df$bind_rows("sample") %>%
    transmute(sample = sample,
              chr = seqnames,
              wgs = ploidy)

names(dimnames(ploidy)) = names(dimnames(wgs)) = c("chr", "sample")
df_inf = reshape2::melt(ploidy) %>%
    mutate(rna = value,
           chr = as.character(chr)) %>%
    select(-value)
dfs = inner_join(df_inf, wgs, by=c("chr", "sample")) %>%
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
             mod %>% broom::glance() %>% pull(r.squared)),
         subtitle = "assuming total DNA=const, diff euploid lib as reference")

pdf("ploidy_eT.pdf")
print(p)
dev.off()

save(ploidy, aneuploidy, wgs, dfs, file="ploidy_eT.RData")
