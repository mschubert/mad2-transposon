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

# ploidy from scWGS -> data.frame with: seqnames, ploidy, sample
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

#' Plot correlation between measured and predicted aneuploidy
#'
#' @param ploidy_df  Data.frame with fields: sample, seqnames, ploidy, ref_ploidy
#' @param ...        Additional parameters passed to labs(...)
#' @return           A ggplot2 object
plot_comparison = function(ploidy_df, ...) {
    mod = lm(ploidy ~ ref_ploidy, data=ploidy_df)
    ploidy_df %>%
        mutate(label=paste(sample, seqnames, sep=":")) %>%
        ggplot(aes(x=ref_ploidy, y=ploidy, label=label, color=sample)) +
        geom_text() +
        geom_smooth(method="lm", color="black") +
        labs(x = "Ploidy from single-cell WGS data + AneuFinder",
             y = "Ploidy inferred from RNA-seq",
             title = sprintf("RNA-seq for ploidy inference (p=%.2g, r^2=%.2f)",
                 mod %>% broom::tidy() %>% filter(term == "ref_ploidy") %>% pull(p.value),
                 mod %>% broom::glance() %>% pull(r.squared)), ...)
}

# control for individual ploidy predictions w/o input sample
ctl = aneuploidy_separate %>%
    filter(sample %in% c("401t", "419t", "413s"),
           sample != ref) %>%
    select(sample, seqnames, ploidy) %>%
    inner_join(kchr %>% select(sample=ref, seqnames, ref_ploidy), by=c("sample", "seqnames"))

# control for min deviation ploidy prediction w/o input sample
ctl2 = ctl %>%
    group_by(seqnames, sample) %>%
    summarize(ploidy = mean(ploidy),
              ref_ploidy = unique(ref_ploidy))

ctl3 = ctl2 %>%
    group_by(sample) %>%
    mutate(ploidy = 2 + ploidy - median(ploidy)) %>%
    ungroup()

pdf("ploidy_cancer_only.pdf")
print(plot_comparison(ctl, subtitle="predicted from one reference sample"))
print(plot_comparison(ctl2, subtitle="predicted from 2 ref samples"))
print(plot_comparison(ctl3, subtitle="predicted from 2 ref samples + median=2"))
dev.off()

#save(ploidy, aneuploidy, wgs, file="ploidy_from_rnaseq.RData")
