library(ggplo2)
library(dplyr)
theme_set(cowplot::theme_cowplot())
io = import('io')
gdsc = import('data/gdsc')
enr = import('tools/enrichr')
idmap = import('process/idmap')
plt = import('plot')
sys = import('sys')

args = sys$cmd$parse(
    opt('p', 'plotfile', 'pdf', 'mile.pdf')
)

clean_reg = function(x) x[!is.na(x) & !duplicated(x)]

dset = io$load("../../data/arrayexpress/E-GEOD-13159.RData")
expr = Biobase::exprs(dset)
rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
meta = Biobase::pData(dset) %>%
    transmute(sample = Array.Data.File,
              type = FactorValue..LEUKEMIA.CLASS.)

# find human blood ets, erg binding in enrichr
tfs = enr$genes("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
mi = io$load("../../data/networks/E-GEOD-13159.RData") %>%
    select(Target, Regulator) %>%
    unstack()
regs = c(tfs[c('ETS1_ENCODE', 'ERG_CHEA')], mi[c('ETS1', 'ERG')])
names(regs) = c("ETS1_chip", "ERG_chip", "ETS1_mi", "ERG_mi")
scores = GSVA::gsva(expr, regs)
idx = cbind(meta, t(scores),
            ERG_expr=expr["ERG",], ETS1_expr=expr["ETS1",])

pdf(args$plotfile, 10, 6)
ggplot(idx, aes(x=ERG_expr, y=ETS1_expr)) +
    geom_point(aes(size=0.5, fill=type), alpha=0.5, shape=21) + # size=aneuploidy
    ggtitle("TF expression")

plt$venn(regs)

ggplot(idx, aes(x=ERG_chip, y=ETS1_chip)) +
    geom_point(aes(size=0.5, fill=type), alpha=0.5, shape=21) + #s=anp
    ggtitle("ChIP regulon enrichment")

ggplot(idx, aes(x=ERG_mi, y=ETS1_mi)) +
    geom_point(aes(size=0.5, fill=type), alpha=0.5, shape=21) + #s=anp
    ggtitle("MI regulon enrichment")

dev.off()
