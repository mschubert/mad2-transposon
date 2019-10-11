library(ggplot2)
library(dplyr)
theme_set(cowplot::theme_cowplot())
io = import('io')
gdsc = import('data/gdsc')
enr = import('tools/enrichr')
plt = import('plot')
sys = import('sys')

args = sys$cmd$parse(
    opt('p', 'plotfile', 'pdf', 'gdsc.pdf')
)

keep = gdsc$tissues(c("PRAD"))
expr = gdsc$basal_expression()
aneup = gdsc$aneuploidy()
narray::intersect(keep, expr, aneup$cosmic, along=-1)

# find human blood ets, erg binding in enrichr
idx = cbind(aneup, tissue=unname(keep),
            ERG_expr=expr["ERG",], ETS1_expr=expr["ETS1",])

pdf(args$plotfile)
ggplot(idx, aes(x=ERG_expr, y=ETS1_expr)) +
    geom_hline(yintercept=3.5, linetype="dashed") +
    geom_vline(xintercept=3.5, linetype="dashed") +
    geom_point(aes(size=aneuploidy, fill=tissue), shape=21) +
    geom_text_repel(aes(label=cell_line), size=2) +
    ggtitle("TF expression")
dev.off()
