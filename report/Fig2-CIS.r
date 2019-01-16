library(dplyr)
library(tidygraph)
library(cowplot)
library(patchwork)
io = import('io')

aneup = io$load("../ploidy_compare/analysis_set.RData") %>%
    select(sample, type, aneuploidy)
ext = io$load("../cis_analysis/ext_gene.RData")
bionet = io$load("../cis_analysis/bionet.RData") %>%
    pull(name)

rna_ins = io$read_table("../data/rnaseq_imfusion/insertions.txt", header=TRUE) %>%
    select(sample, external_gene_name=gene_name) %>%
    filter(toupper(external_gene_name) %in% bionet) %>%
    distinct() %>%
    mutate(rna_ins = 1)

cis = io$load("../cis_analysis/poisson_gene.RData")
cis_result = cis$result %>%
    ungroup() %>% # TODO: don't saved grouped
    filter(toupper(external_gene_name) %in% bionet,
           adj.p < 1e-3)
cis_samples = cis$samples %>%
    filter(external_gene_name %in% cis_result$external_gene_name) %>%
    tidyr::complete(sample, external_gene_name, fill=list(n_ins=0, reads=0)) %>%
    group_by(sample) %>%
    mutate(total_ins = sum(n_ins),
           total_reads = sum(reads),
           gene_read_frac = reads / total_reads) %>%
    ungroup() %>%
    left_join(rna_ins) %>%
    mutate(rna_ins = ifelse(is.na(rna_ins), 0, 1),
        ins_type = case_when(
            n_ins > 0 & rna_ins > 0 ~ "both",
            n_ins > 0 & rna_ins == 0 ~ "DNA",
            n_ins == 0 & rna_ins > 0 ~ "RNA"
        ),
        ins_type = factor(ins_type, levels=c("DNA", "RNA", "both"))
    ) %>%
    inner_join(aneup) %>%
    inner_join(ext$aneuploidy %>% select(external_gene_name, aneup_stat=statistic)) %>%
    mutate(sample = forcats::fct_reorder(sample, aneuploidy),
           external_gene_name = forcats::fct_reorder(external_gene_name, aneup_stat))
genelvl = levels(cis_samples$external_gene_name)
smplvl = levels(cis_samples$sample)
cis_stats = cis_result %>%
    mutate(external_gene_name = factor(external_gene_name, levels=genelvl)) %>%
    filter(!is.na(external_gene_name))

types = bind_rows(ext) %>%
    filter(subset %in% c("Myeloid", "T-cell", "Other")) %>%
    filter(toupper(external_gene_name) %in% bionet) %>%
    arrange(p.value) %>%
    group_by(subset) %>%
    top_n(4, -p.value)

p1 = ggplot(cis_samples, aes(x=sample, y=external_gene_name)) +
    geom_tile(aes(fill=ins_type, alpha=gene_read_frac)) +
    scale_fill_manual(values=c("maroon4", "navy", "springgreen4")) +
    theme(axis.text.x = element_text(angle=90, vjust=0.5))

p11 = mutate(aneup, sample=factor(sample, smplvl)) %>%
    filter(!is.na(sample)) %>%
    ggplot(aes(x=sample, y=aneuploidy, fill=type)) +
    geom_bar(stat="identity") +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank())

p12 = ggplot(cis_stats, aes(x=external_gene_name, y=-log10(adj.p))) +
    geom_bar(stat="identity") +
    scale_x_discrete(position="top") +
    coord_flip() +
    theme(axis.title.y = element_blank())

p11 + plot_spacer() + p1 + p12 + plot_layout(widths=c(5,1), heights=c(1,5))

p2 = ggplot(types, aes(x=forcats::fct_reorder(external_gene_name, statistic), y=statistic)) +
    geom_hline(yintercept=0, color="grey", linetype="dashed") +
    geom_bar(aes(fill=subset), stat="identity") +
    geom_text(aes(label=external_gene_name, y=statistic/2)) +
    coord_flip() +
    facet_wrap(~ subset, scales="free_y") +
    theme(axis.text.y = element_blank(),
          axis.line = element_blank(),
          axis.ticks.y = element_blank()) +
    guides(fill = FALSE) +
    labs(y = "Wald statistic",
         x = "Gene")


