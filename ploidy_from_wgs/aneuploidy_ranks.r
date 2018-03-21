library(dplyr)
library(cowplot)
library(patchwork)
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'tsv with aneuploidy scores', 'copy_segments.tsv'),
    opt('p', 'plotfile', 'pdf to save plot to', '/dev/null'))

aneup = readr::read_tsv(args$infile) %>%
    mutate(Sample = forcats::fct_reorder(Sample, aneuploidy))

dens = ggplot(aneup, aes(x=aneuploidy)) +
    geom_density(alpha=0.5, fill="red") +
    theme(axis.title.x=element_blank(),
          axis.text=element_blank(),
          axis.line=element_blank(),
          axis.ticks=element_blank())

samples = ggplot(aneup, aes(x=aneuploidy, y=Sample)) +
#    geom_vline(xintercept=0.05, color="blue", linetype="dashed") +
    geom_segment(aes(xend=aneuploidy, y=Sample, yend=Sample), x=0, color="lightgrey") +
    geom_point() +
    theme(axis.text.y = element_text(size=8))

pdf(4, 10, file=args$plotfile)
print(dens + samples + plot_layout(ncol=1, heights=c(1,12)))
dev.off()
