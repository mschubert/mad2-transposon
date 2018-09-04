library(cowplot)
library(dplyr)
io = import('io')
gdsc = import('data/gdsc')
enr = import('tools/enrichr')
plt = import('plot')

keep = gdsc$tissues(c("LAML", "DLBC", "LCML", "ALL"))
expr = gdsc$basal_expression()
aneup = gdsc$aneuploidy()
narray::intersect(keep, expr, aneup$cosmic, along=-1)

# find human blood ets, erg binding in enrichr
tfs = enr$genes("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
mi = io$load("../../data/networks/E-GEOD-13159.RData") %>%
    select(Target, Regulator) %>%
    unstack()
regs = c(tfs[c('ETS1_ENCODE', 'ERG_CHEA')], mi[c('ETS1', 'ERG')])
names(regs) = c("ETS1_chip", "ERG_chip", "ETS1_mi", "ERG_mi")
scores = GSVA::gsva(expr, regs)
idx = cbind(aneup, tissue=unname(keep), t(scores),
            ERG_expr=expr["ERG",], ETS1_expr=expr["ETS1",])

pdf("gdsc.pdf")
ggplot(idx, aes(x=ERG_expr, y=ETS1_expr)) +
    geom_hline(yintercept=3.5, linetype="dashed") +
    geom_vline(xintercept=3.5, linetype="dashed") +
    geom_point(aes(size=aneuploidy, fill=tissue), shape=21) +
    geom_text_repel(aes(label=cell_line), size=2) +
    ggtitle("TF expression")

plt$venn(regs)

ggplot(idx, aes(x=ERG_chip, y=ETS1_chip)) +
    geom_point(aes(size=aneuploidy, fill=tissue), shape=21) +
    geom_text_repel(aes(label=cell_line), size=2) +
    ggtitle("ChIP regulon enrichment")

idx2 = idx %>%
    filter(ERG_expr > 3.5 & ETS1_expr > 3.5, tissue == "ALL")
idx2$aneup_class = ifelse(idx2$aneuploidy > 0.5, "aneuploid", "euploid")
p = lm(ETS1_chip ~ aneup_class:ERG_chip, data=idx2) %>%
    broom::tidy() %>% pull(p.value)
ggplot(idx2, aes(x=ERG_chip, y=ETS1_chip, fill=aneup_class)) +
    geom_point(aes(size=aneuploidy), shape=21) +
    geom_smooth(aes(color=aneup_class), method="lm", se=FALSE) +
    geom_text_repel(aes(label=cell_line), size=2) +
    ggtitle(sprintf("ALL regulon enrichment, TF expressed (p eup %.2f, anp %.2f)", p[3], p[2]))

ggplot(idx, aes(x=ERG_mi, y=ETS1_mi)) +
    geom_point(aes(size=aneuploidy, fill=tissue), shape=21) +
    geom_text_repel(aes(label=cell_line), size=2) +
    ggtitle("MI regulon enrichment")

dev.off()
