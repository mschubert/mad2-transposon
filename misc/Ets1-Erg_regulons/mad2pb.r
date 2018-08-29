library(cowplot)
library(dplyr)
io = import('io')
enr = import('tools/enrichr')
idmap = import('process/idmap')

dset = io$load("../../data/rnaseq/assemble.RData")
expr = dset$expr
rownames(expr) = dset$genes
meta = io$load("../../ploidy_compare/analysis_set.RData")
narray::intersect(expr, meta$sample, along=-1)

# find human blood ets, erg binding in enrichr
tfs = stack(enr$genes("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X"))
tfs$values = unname(idmap$orthologue(tfs$values, from="external_gene_name", to="mgi_symbol"))
tfs = unstack(tfs)
scores = GSVA::gsva(expr, tfs)

idx = cbind(meta %>% select(sample, aneuploidy, type),
            Erg_expr=expr["Erg",], Ets1_expr=expr["Ets1",],
            Erg_reg=scores["ERG_CHEA",], Ets1_reg=scores["ETS1_ENCODE",])

pdf("mad2pb.pdf")
ggplot(idx, aes(x=Erg_expr, y=Ets1_expr)) +
    geom_point(aes(size=aneuploidy, fill=type), shape=21) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    ggtitle("TF expression")

ggplot(idx, aes(x=Erg_reg, y=Ets1_reg)) +
    geom_point(aes(size=aneuploidy, fill=type), shape=21) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    ggtitle("Regulon enrichment")

dev.off()
