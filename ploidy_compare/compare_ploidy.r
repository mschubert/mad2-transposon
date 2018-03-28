library(dplyr)
library(cowplot)
library(patchwork)
io = import('io')

#args = sys$cmd$parse(
#    opt('i', 'infile', 'RData aneuploidy scores', 'copy_segments.RData'),
#    opt('p', 'plotfile', 'pdf to save plot to', '/dev/null'))

dna = io$load("../ploidy_from_wgs/copy_segments.RData")$aneups %>%
    dplyr::rename(sample = Sample) %>%
    mutate(sample = sub("(-high)|(-low)", "", sample))
rna = io$load("../ploidy_from_rnaseq/panel_ploidy.RData")$aneup %>%
    mutate(sample = toupper(gsub("[^0-9stST]+", "", sample)))

aneups = list(WGS = dna, `RNA-seq` = rna) %>%
    dplyr::bind_rows(.id="type") %>%
    mutate(sample = forcats::fct_reorder(sample, aneuploidy))

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
    geom_point(aes(color=type)) +
    theme(axis.text.y = element_text(size=8))

pdf(7, 9, file="compare_ploidy.pdf")
print(dens + samples + plot_layout(ncol=1, heights=c(1,12)))
dev.off()
