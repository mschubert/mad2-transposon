library(cowplot)
library(dplyr)
io = import('io')
ccle = import('data/ccle')
plt = import('plot')
sys = import('sys')

args = sys$cmd$parse(
    opt('p', 'plotfile', 'pdf', 'Ets1+Erg_ccle.pdf')
)

meta = as.data.frame(ccle$index)
keep = ccle$tissues(c("haematopoietic_and_lymphoid_tissue"))
expr = ccle$basal_expression()
#aneup = ccle$aneuploidy()
narray::intersect(keep, expr, meta$`Expression arrays`, along=-1)

# find human blood ets, erg binding in enrichr
idx = data.frame(cell_line = meta$`Cell line primary name`,
                 tissue = meta$`Hist Subtype1`,
                 ERG_expr = expr["ERG",],
                 ETS1_expr = expr["ETS1",])

pdf(args$plotfile)
ggplot(idx, aes(x=ERG_expr, y=ETS1_expr)) +
    geom_hline(yintercept=3.5, linetype="dashed") +
    geom_vline(xintercept=3.5, linetype="dashed") +
    geom_point(aes(size=1, fill=tissue), shape=21) +
    geom_text_repel(aes(label=cell_line), size=2) +
    ggtitle("TF expression")
dev.off()
