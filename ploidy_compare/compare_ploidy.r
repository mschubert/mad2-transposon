library(dplyr)
library(cowplot)
library(patchwork)
io = import('io')

#args = sys$cmd$parse(
#    opt('i', 'infile', 'RData aneuploidy scores', 'copy_segments.RData'),
#    opt('p', 'plotfile', 'pdf to save plot to', '/dev/null'))

dna = io$load("../ploidy_from_wgs/copy_segments.RData")
rna = io$load("../ploidy_from_rnaseq/panel_ploidy.RData")

aneups = list(
        WGS = dna$aneups %>% rename(sample = Sample),
        `RNA-seq` = rna$aneups %>%
            mutate(sample = toupper(gsub("[^0-9ST]", "", sample)))) %>%
    dplyr::bind_rows("type") %>%
    mutate(sample = forcats::fct_reorder(sample, aneuploidy))

dens = ggplot(aneup, aes(x=aneuploidy)) +
    geom_density(aes(fill="red"), alpha=0.5) +
    theme(axis.title.x=element_blank(),
          axis.text=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank())

samples = ggplot(aneup, aes(x=aneuploidy, y=sample, color=type)) +
#    geom_vline(xintercept=0.05, color="blue", linetype="dashed") +
    geom_segment(aes(xend=aneuploidy, y=sample, yend=sample), x=0, color="lightgrey") +
    geom_point() +
    theme(axis.text.y = element_text(size=8))

pdf(4, 10, file=args$plotfile)
print(dens + samples + plot_layout(ncol=1, heights=c(1,12)))
dev.off()
