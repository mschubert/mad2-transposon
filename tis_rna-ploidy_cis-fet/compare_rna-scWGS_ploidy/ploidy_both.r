library(cowplot)
library(ggrepel)
io = import('io')
df = import('data_frame')
idmap = import('process/idmap')
aneufinder = import('tools/aneufinder')

# reads per chromosome for T-ALLs
counts1 = io$load('../../data/rnaseq/assemble.RData')$counts
expr2 = readr::read_tsv("../../data/rnaseq/T-ALL_read_count.txt")
counts2 = data.matrix(expr2[,c("eT_p0","eT_p2")])
counts = narray::stack(list(counts1, counts2), along=2) %>% na.omit()

chrs = idmap$gene(rownames(counts), to="chromosome_name", dset="mmusculus_gene_ensembl")
reads = counts[!is.na(chrs),]
chrs = chrs[!is.na(chrs)]
reads = narray::map(reads, along=1, sum, subsets=chrs)[as.character(1:19),]
#reads = reads / narray::rrep(colSums(reads), nrow(reads))
names(dimnames(reads)) = c("seqnames", "sample")
read_df = reshape2::melt(reads, value.name = "reads") %>%
    mutate(seqnames = as.character(seqnames),
           sample = as.character(sample))

# ploidy from scWGS -> data.frame with: seqnames, ploidy, sample
wgs = c("T401", "T419", "S413") %>%
    paste0("../../../dosage/data/singlecell_wgs/T-ALL/", ., ".RData") %>%
    io$load() %>%
    lapply(aneufinder$consensus_ploidy) %>%
    setNames(c("401t", "419t", "413s")) %>%
    df$bind_rows("sample")

add = expand.grid(seqnames=unique(wgs$seqnames), ploidy=2, sample=c("eT_p0", "eT_p2"))
wgs = dplyr::bind_rows(wgs, add)

eff_libsize = wgs %>%
    dplyr::left_join(read_df, by=c("sample", "seqnames")) %>%
    na.omit() %>% # 'X' missing
    group_by(sample) %>%
    summarize(libsize_per_somy = sum(reads / ploidy))

eff_reads = dplyr::inner_join(read_df, eff_libsize, by="sample") %>%
    mutate(eff_reads = reads / libsize_per_somy) %>%
    select(-libsize_per_somy, -reads)

# calculate k^chr (cf. chr_ngenes_nreads_cor.r)
kchr = wgs %>%
    dplyr::left_join(eff_reads, by=c("sample", "seqnames")) %>%
    mutate(kchr = ploidy / eff_reads) %>%
    transmute(seqnames=seqnames, ref=sample, ref_ploidy=ploidy, kchr=kchr) %>%
    arrange(seqnames, ref)

# aneuploidy scores for all samples
aneuploidy_separate = read_df %>%
    left_join(kchr, by="seqnames") %>%
    mutate(ploidy = reads * kchr) %>%
    group_by(sample) %>%
    mutate(ploidy = 2 * ploidy / median(ploidy)) %>%
    ungroup()

#aneuploidy_mean = aneuploidy_separate %>%
#    group_by(sample, seqnames) %>%
#    summarize(frac_reads = unique(frac_reads),
#              ploidy = mean(ploidy))

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
    mutate(ploidy = ploidy - 2) %>%
    group_by(seqnames, sample) %>%
    summarize(ploidy = mean(ploidy),
              ref_ploidy = unique(ref_ploidy))

ctl3 = ctl2 %>%
    group_by(sample) %>%
    mutate(ploidy = 2 + ploidy - median(ploidy)) %>%
    ungroup()

pdf("ploidy_both.pdf")
print(plot_comparison(ctl, subtitle="predicted from one reference sample"))
print(plot_comparison(ctl2, subtitle="predicted from 2 ref + 2 eT samples"))
print(plot_comparison(ctl3, subtitle="predicted from 2 ref + 2 eT samples + median=2"))
dev.off()

save(ctl, ctl2, ctl3, file="ploidy_both.RData")
