library(cowplot)
library(dplyr)
gdsc = import('data/gdsc')
enr = import('tools/enrichr')

keep = gdsc$tissues(c("LAML", "DLBC", "LCML", "ALL"))
expr = gdsc$basal_expression()
aneup = gdsc$aneuploidy()
narray::intersect(keep, expr, aneup$cosmic, along=-1)

# find human blood ets, erg binding in enrichr
tfs = enr$genes("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
scores = GSVA::gsva(expr, tfs)

idx = cbind(aneup, tissue=unname(keep),
            ERG_expr=expr["ERG",], ETS1_expr=expr["ETS1",],
            ERG_reg=scores["ERG_CHEA",], ETS1_reg=scores["ETS1_ENCODE",])

pdf("gdsc.pdf")
ggplot(idx, aes(x=ERG_expr, y=ETS1_expr)) +
    geom_hline(yintercept=3.5, linetype="dashed") +
    geom_vline(xintercept=3.5, linetype="dashed") +
    geom_point(aes(size=aneuploidy, fill=tissue), shape=21) +
    geom_text_repel(aes(label=cell_line), size=2) +
    ggtitle("TF expression")

ggplot(idx, aes(x=ERG_reg, y=ETS1_reg)) +
    geom_point(aes(size=aneuploidy, fill=tissue), shape=21) +
    geom_text_repel(aes(label=cell_line), size=2) +
    ggtitle("Regulon enrichment")

idx2 = idx %>%
    filter(ERG_expr > 3.5 & ETS1_expr > 3.5, tissue == "ALL")
idx2$aneup_class = ifelse(idx2$aneuploidy > 0.5, "aneuploid", "euploid")
p = lm(ETS1_reg ~ aneup_class:ERG_reg, data=idx2) %>%
    broom::tidy() %>% pull(p.value)
ggplot(idx2, aes(x=ERG_reg, y=ETS1_reg, fill=aneup_class)) +
    geom_point(aes(size=aneuploidy), shape=21) +
    geom_smooth(aes(color=aneup_class), method="lm", se=FALSE) +
    geom_text_repel(aes(label=cell_line), size=2) +
    ggtitle(sprintf("ALL regulon enrichment, TF expressed (p eup %.2f, anp %.2f)", p[3], p[2]))
dev.off()
