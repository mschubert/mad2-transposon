library(dplyr)
st = import('stats')
plt = import('plot')
enr = import('tools/enrichr')
idmap = import('process/idmap')

# response signatures: Stat1 ChIP (enrichr), Ifn response
expr = readRDS("../expr_diff/eset_Mad2PB.rds")
#genes = c("Stat1", "Pias1", "Mb21d1", "Ifng", "Ifngr1", "Trp53", "Wrap53")
genes = c("Stat1", "Pias1", "Ifng", "Trp53")
ee = readRDS("../expr_diff/eset_Mad2PB.rds")$vs
stat1 = ee["Stat1",]
rest = ee[-which(rownames(ee) == "Stat1"),]
cc = data.frame(
    gene = rownames(rest),
    cor = narray::map(rest, along=2, function(x) cor(x, stat1))
) %>% arrange(-cor)

encc = readRDS("../data/genesets/mouse/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.rds")
chea = readRDS("../data/genesets/mouse/ChEA_2016.rds")
go = readRDS("../data/genesets/mouse/GO_Biological_Process_2018.rds")
hm = readRDS("../data/genesets/mouse/CH.HALLMARK.rds")
mile = readRDS("../data/genesets/mouse/MILE_regulons.rds")

sets=c(chea[c("STAT1_17558387_ChIP-Seq_HELA_Human",
              "STAT3_23295773_ChIP-Seq_U87_Human")],
       hm["HALLMARK_INTERFERON_GAMMA_RESPONSE"])
sets = list(
    IFNg_Stat1 = intersect(sets[[1]], sets[[3]]),
    IFNg_Stat3 = intersect(sets[[2]], sets[[3]]),
    IFNg_other = setdiff(sets[[3]], c(sets[[1]], sets[[2]])),
         STAT1_cor = setdiff(intersect(
            head(cc$gene, 1000),
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]]), sets[[3]])
)

#sets = c(go["regulation of complement activation (GO:0030449)"],
#         chea["STAT1_20625510_ChIP-Seq_HELA_Human"],
#         chea["STAT1_17558387_ChIP-Seq_HELA_Human"],
#         encc["STAT3_CHEA"],
#         encc["STAT3_ENCODE"],
#         hm["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
#         STAT1_complement = list(intersect(
#            go[["regulation of complement activation (GO:0030449)"]],
#            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
#         STAT1_mhc = list(intersect(
#            go[["antigen receptor-mediated signaling pathway (GO:0050851)"]],
#            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
#         STAT1_apop = list(intersect(
#            go[["apoptotic process (GO:0006915)"]],
#            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
#         STAT1_cor = list(intersect(
#            head(cc$gene, 1000),
#            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]]))
#)
#sets$HALLMARK_INTERFERON_GAMMA_RESPONSE = setdiff(sets$HALLMARK_INTERFERON_GAMMA_RESPONSE,
#                                                  chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])
## + MHC, IFN(g), apop(?) -> do GO enrichment, then choose (or: clust?)
gsva = GSVA::gsva(expr$vs, sets)
scores = rbind(expr$vs[genes,], gsva)

# expression changes with insertions (in different subtypes)
cis = readRDS("../cis_analysis/poisson.rds")$samples %>%
    filter(external_gene_name %in% genes, sample %in% colnames(expr$vs)) %>%
    narray::construct(n_ins ~ sample + external_gene_name, data=., fill=0) %>%
    apply(1, function(x) paste(names(x)[x != 0], collapse="+"))
annot = as.data.frame(SummarizedExperiment::colData(expr$eset))
annot$ins = "none"
annot$ins[match(names(cis), annot$sample)] = cis
annot$ins = relevel(factor(annot$ins), "none")

both = cbind(annot, t(scores)) %>%
    filter(type != "unknown")

cors = both %>% 
    tidyr::gather("subs", "expr", IFNg_Stat1:IFNg_other)

cors2 = cors %>%
#    filter(subs != "STAT1_mhc") %>%
    mutate(subs = factor(subs))
#levels(cors2$subs) = c("Apoptotic process\n(GO:0006915)",
#                       "Regulation of complement\nactivation (GO:0030449)")
p2 = ggplot(cors2, aes(x=STAT1_cor, y=expr, color=type)) +
    geom_point(aes(size=aneuploidy), alpha=0.8) +
    geom_smooth(method='lm', color="black", se=FALSE) +
    geom_smooth(method='lm', aes(color=type), se=FALSE) +
    coord_fixed() +
    facet_wrap(~ subs) +
    guides(color = guide_legend(title="Cancer type"),
           size = guide_legend(title="Aneuploidy")) +
    labs(x = "Stat1 activity (inferred)",
         y = "GO signature score")









# add cor plots for stat1 act <> stat1 subsets (complement, mhc, apop)
# + <> aneup
# pt size: stat1 expr?
#    tidyr::gather("subs", "expr", aneuploidy, STAT1_complement:STAT1_apop)

stats = split(cors, cors$subs) %>%
    lapply(function(x) broom::tidy(lm(expr ~ STAT1_cor, data=x))) %>%
    dplyr::bind_rows(.id="with") %>%
    filter(term == "STAT1_cor") %>%
    select(-term)
stats2 = cors %>%
    filter(type == "Other") %>%
    split(.$subs) %>%
    lapply(function(x) broom::tidy(lm(expr ~ STAT1_cor, data=x))) %>%
    dplyr::bind_rows(.id="with") %>%
    filter(term == "STAT1_cor") %>%
    select(-term)
stats3 = broom::tidy(lm(aneuploidy ~ STAT1_cor, data=cors)) %>%
    filter(term == "STAT1_cor") %>%
    select(-term)
stats4 = broom::tidy(lm(aneuploidy ~ STAT1_cor, data=filter(cors, type=="Other"))) %>%
    filter(term == "STAT1_cor") %>%
    select(-term)

p1 = ggplot(cors, aes(x=STAT1_cor, y=aneuploidy)) +
    geom_vline(xintercept=0, linetype="dashed") +
    geom_boxplot(aes(group=Stat1_act), outlier.shape=NA, color="grey", varwidth=TRUE) +
    geom_boxplot(aes(group=Stat1_act, color=type), outlier.shape=NA, color="grey", varwidth=TRUE) +
    geom_point(aes(color=type), size=3) +
    guides(color = guide_legend(title="Cancer type")) +
    labs(x = "Stat1 activity (inferred)",
         y = "Aneuploidy")
cors2 = cors %>%
#    filter(subs != "STAT1_mhc") %>%
    mutate(subs = factor(subs))
#levels(cors2$subs) = c("Apoptotic process\n(GO:0006915)",
#                       "Regulation of complement\nactivation (GO:0030449)")
p2 = ggplot(cors2, aes(x=STAT1_cor, y=expr, color=type)) +
    geom_point(aes(size=aneuploidy), alpha=0.8) +
    geom_smooth(method='lm', color="black", se=FALSE) +
    geom_smooth(method='lm', aes(color=type), se=FALSE) +
    coord_fixed() +
    facet_wrap(~ subs) +
    guides(color = guide_legend(title="Cancer type"),
           size = guide_legend(title="Aneuploidy")) +
    labs(x = "Stat1 activity (inferred)",
         y = "GO signature score")

library(gridExtra)
pdf("stat1_compare_stat346.pdf", 10, 8)
print(p1)
print(p2)
grid.arrange(grid::textGrob("Aneuploidy Stat1 low vs high, all types"), tableGrob(stats3))
grid.arrange(grid::textGrob("Aneuploidy Stat1 low vs high, B-like only"), tableGrob(stats4))
grid.arrange(grid::textGrob("Stat1+GO line fit, all types"), tableGrob(stats))
grid.arrange(grid::textGrob("Stat1+GO line fit, B-like only"), tableGrob(stats2))
dev.off()
