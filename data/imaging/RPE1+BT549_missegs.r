library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'xlsx', '20190819 analysis movie 20190815.xlsx'),
    opt('p', 'plotfile', 'pdf', 'RPE1+BT549_missegs.pdf'))

cols = rev(c('#7fc97f', '#beaed4', '#fdc086', '#ffff99', '#386cb0', '#f0027f',
             '#bf5b17', '#666666'))

dset = list(
    rpe1 = readxl::read_xlsx(args$infile, "RPE-1 H2B GFP"),
    bt549 = readxl::read_xlsx(args$infile, "BT549 H2B-GFP")
) %>% bind_rows(.id="cell_line") %>%
    mutate(mitotic_time = `mitotic exit` - `mitotic entry`,
           mitotic_exit = `real time (min)` + mitotic_time) %>%
    select(-...7, -`real time`) %>%
    na.omit()
colnames(dset) = make.names(colnames(dset))
dset$type.of.mitosis = relevel(factor(dset$type.of.mitosis), "normal")

# plot read-like stacks of mitotic timing
library(GenomicRanges)
library(ggbio)
gr = makeGRangesFromDataFrame(dset, keep.extra.columns=TRUE, seqnames.field="cell_line",
    start.field="real.time..min.", end.field="mitotic_exit")
#ggbio::plotStackedOverview(gr)

p = autoplot(gr, aes(fill=type.of.mitosis)) +
    facet_wrap(~ seqnames, ncol=1) +
    scale_fill_manual(values=cols) +
    scale_x_continuous("hours", breaks=c(0:10)*600, labels=c(0:10)*600/60) +
    geom_vline(xintercept=c(24,48,72)*60, linetype="dashed") +
    ggtitle("mitoses rev 500 nM")

p2 = ggplot(dset, aes(x=real.time..min., fill=type.of.mitosis)) +
    facet_wrap(~ cell_line, ncol=1) +
    ggridges::geom_density_ridges(aes(y=type.of.mitosis, height=..count..),
                                  alpha=0.5, stat="binline", bins=15, scale=4) +
#    geom_density(alpha=0.3) +
    scale_fill_manual(values=setNames(cols, levels(dset$type.of.mitosis))) +
    scale_x_continuous("hours", breaks=c(0:10)*600, labels=c(0:10)*600/60, limits=c(7,78)*60) +
    geom_vline(xintercept=c(24,48,72)*60, linetype="dashed") +
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank())

pdf(args$plotfile, 15, 3)
print(p)
print(p2)
dev.off()
