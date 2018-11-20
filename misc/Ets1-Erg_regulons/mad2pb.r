library(cowplot)
library(dplyr)
io = import('io')
gdsc = import('data/gdsc')
enr = import('tools/enrichr')
idmap = import('process/idmap')
plt = import('plot')
sys = import('sys')

args = sys$cmd$parse(
    opt('p', 'plotfile', 'pdf', 'mad2pb.pdf')
)

clean_reg = function(x) x[!is.na(x) & !duplicated(x)]

dset = io$load("../../data/rnaseq/assemble.RData")
expr = dset$expr
rownames(expr) = dset$genes
meta = io$load("../../ploidy_compare/analysis_set.RData") %>%
    select(sample, aneuploidy, type)
narray::intersect(expr, meta$sample, along=-1)

# find human blood ets, erg binding in enrichr
tfs = enr$genes("ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X")
mi = io$load("../../data/networks/E-GEOD-13159.RData") %>%
    select(Target, Regulator) %>%
    unstack()
regs = c(tfs[c('ETS1_ENCODE', 'ERG_CHEA')], mi[c('ETS1', 'ERG')]) %>%
    idmap$orthologue(from="external_gene_name", to="mgi_symbol") %>%
    lapply(clean_reg)
names(regs) = c("ETS1_chip", "ERG_chip", "ETS1_mi", "ERG_mi")
scores = GSVA::gsva(expr, regs)
idx = cbind(meta, t(scores),
            Erg_expr=expr["Erg",], Ets1_expr=expr["Ets1",],
            Etv6_expr=expr["Etv6",], Elf1_expr=expr["Elf1",])

pdf(args$plotfile)
ggplot(idx, aes(x=Erg_expr, y=Ets1_expr)) +
    geom_point(aes(size=aneuploidy, fill=type), shape=21) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    ggtitle("TF expression")

ggplot(idx, aes(x=Etv6_expr, y=Elf1_expr)) +
    geom_point(aes(size=aneuploidy, fill=type), shape=21) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    ggtitle("TF expression")

ggplot(idx, aes(x=Etv6_expr, y=Erg_expr)) +
    geom_point(aes(size=aneuploidy, fill=type), shape=21) +
    ggrepel::geom_text_repel(aes(label=sample), size=2) +
    ggtitle("TF expression")

plt$venn(regs)

ggplot(idx, aes(x=ERG_chip, y=ETS1_chip)) +
    geom_point(aes(size=aneuploidy, fill=type), shape=21) +
    geom_text_repel(aes(label=sample), size=2) +
    ggtitle("ChIP regulon enrichment")

ggplot(idx, aes(x=ERG_mi, y=ETS1_mi)) +
    geom_point(aes(size=aneuploidy, fill=type), shape=21) +
    geom_text_repel(aes(label=sample), size=2) +
    ggtitle("MI regulon enrichment")

dev.off()
