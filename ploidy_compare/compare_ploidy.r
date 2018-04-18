library(dplyr)
library(cowplot)
library(patchwork)
io = import('io')
seq = import('seq')

#args = sys$cmd$parse(
#    opt('i', 'infile', 'RData aneuploidy scores', 'copy_segments.RData'),
#    opt('p', 'plotfile', 'pdf to save plot to', '/dev/null'))

dna = io$load("../ploidy_from_wgs/copy_segments.RData")$segments %>%
    seq$aneuploidy() %>%
    dplyr::rename(sample = Sample) %>%
    mutate(sample = sub("(-high)|(-low)", "", sample),
           sample = paste0(sub("[^ST]+", "", sample), sub("[^0-9]+", "", sample)))
rna = io$load("../ploidy_from_rnaseq/panel_ploidy.RData")$segments %>%
    seq$aneuploidy(sample="sample", ploidy="expr") %>%
    mutate(sample = toupper(gsub("[^0-9stST]+", "", sample)))
eT2 = io$load("../ploidy_from_rnaseq/eT_ploidy.RData")$segments %>%
    seq$aneuploidy(sample="sample", ploidy="expr") %>%
    mutate(sample = toupper(gsub("[^0-9stST]+", "", sample)))
eT = io$load("../tis_rna-ploidy_cis-fet/aneuploidy_mad2.RData") %>%
    transmute(sample=sample, aneuploidy=aneup)

tissues = setNames(c("spleen", "thymus"), c("S","T"))

#aneups = list(WGS = dna, `RNA-seq (eDivisive)` = rna, `RNA-seq (eT)` = eT, `RNA-seq (eT, new)` = eT2) %>%
aneups = list(WGS = dna, `RNA-seq (eT)` = eT2) %>%
    dplyr::bind_rows(.id="type") %>%
    mutate(sample = forcats::fct_reorder(sample, aneuploidy),
           tissue = tissues[sub("Healthy|[0-9]+", "", sample)])

dens = ggplot(aneups, aes(x=aneuploidy)) +
    geom_density(aes(fill=type), alpha=0.5) +
    theme(axis.title.x=element_blank(),
          axis.text=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank()) +
    guides(fill=FALSE)

samples = ggplot(aneups, aes(x=aneuploidy, y=sample)) +
    geom_segment(aes(xend=aneuploidy, yend=sample), x=0,
                 color="lightgrey") +
    geom_point(aes(color=type, shape=tissue)) +
    theme(axis.text.y = element_text(size=8))

paircor = aneups %>%
    group_by(sample, type, tissue) %>%
    summarize(aneuploidy = mean(aneuploidy)) %>%
    tidyr::spread("type", "aneuploidy") %>%
    GGally::ggpairs(columns=3:ncol(.), aes(shape=tissue))

pdf(9, 9, file="compare_ploidy.pdf")
print(dens + samples + plot_layout(ncol=1, heights=c(1,12)))
print(paircor)
dev.off()
