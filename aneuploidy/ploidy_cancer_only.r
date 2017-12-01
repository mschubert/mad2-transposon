library(cowplot)
library(ggrepel)
io = import('io')
df = import('data_frame')
idmap = import('process/idmap')
aneufinder = import('../../aneuploidy/data/singlecell_wgs/aneufinder')

# reads per chromosome for T-ALLs
dset = io$load('../data/rnaseq/assemble.RData')
chrs = idmap$gene(rownames(dset$counts), to="chromosome_name", dset="mmusculus_gene_ensembl")
reads = dset$counts[!is.na(chrs),]
chrs = chrs[!is.na(chrs)]
reads = narray::map(reads, along=1, sum, subsets=chrs)[as.character(1:19),]
reads = reads / narray::rrep(colSums(reads), nrow(reads))
names(dimnames(reads)) = c("seqnames", "sample")
read_df = reshape2::melt(reads, value.name = "frac_reads") %>%
    mutate(seqnames = as.character(seqnames),
           sample = as.character(sample))

# ploidy from scWGS
wgs = c("T401", "T419", "S413") %>%
    paste0("../../aneuploidy/data/singlecell_wgs/T-ALL/", ., ".RData") %>%
    io$load() %>%
    lapply(aneufinder$consensus_ploidy) %>%
    setNames(c("401t", "419t", "413s")) %>%
    df$bind_rows("sample")

# calculate k^chr (cf. chr_ngenes_nreads_cor.r)
kchr = wgs %>%
    dplyr::left_join(read_df, by=c("sample", "seqnames")) %>%
    mutate(kchr = ploidy / frac_reads) %>%
    transmute(seqnames=seqnames, ref=sample, ref_ploidy=ploidy, kchr=kchr)

# aneuploidy scores for all samples
aneuploidy_separate = read_df %>%
    left_join(kchr, by="seqnames") %>%
    mutate(ploidy = frac_reads * kchr)

aneuploidy_mean = aneuploidy_separate %>%
    group_by(sample, seqnames) %>%
    summarize(frac_reads = unique(frac_reads),
              ploidy = mean(ploidy))

#aneuploidy_min = aneuploidy_separate %>%
#    group_by(sample, seqnames) %>%
#    summarize(frac_reads = unique(frac_reads),
#              ploidy = sign(ploidy) * (2 + min(abs(ploidy)-2)))

# control for individual ploidy predictions w/o input sample
ctl = aneuploidy_separate %>%
    filter(sample %in% c("401t", "419t", "413s"),
           sample != ref)

mod = lm(ploidy ~ ref_ploidy, data=ctl)
p_sep = ggplot(ctl, aes(x=ref_ploidy, y=ploidy, label=seqnames, color=ref)) +
    geom_text() +
    geom_smooth(method="lm", color="black") +
    labs(x = "Ploidy from single-cell WGS data + AneuFinder",
         y = "Ploidy inferred from RNA-seq",
         title = sprintf("RNA-seq for ploidy inference (p=%.2g, r^2=%.2f)",
             mod %>% broom::tidy() %>% filter(term == "ref_ploidy") %>% pull(p.value),
             mod %>% broom::glance() %>% pull(r.squared)))

# control for min deviation ploidy prediction w/o input sample

#TODO: calc min deviation from <ctl> df

pdf("ploidy_cancer_only.pdf")
print(p_sep)
#print(p_min)
dev.off()

save(ploidy, aneuploidy, wgs, file="ploidy_from_rnaseq.RData")
