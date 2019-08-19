library(dplyr)
io = import('io')
tcga = import('data/tcga')
idmap = import('process/idmap')
util = import('../../cgas/cor_structure/genenet')

gsva_tcga = function(cohort, genes, sets) {
    expr = tcga$rna_seq(cohort, trans="raw")
    keep = rowMeans(expr, na.rm=TRUE) >= 10
    expr = tcga$rna_seq(cohort, trans="vst")
    rownames(expr) = idmap$gene(rownames(expr), to="external_gene_name")
    expr = expr[keep | rownames(expr) %in% genes,]
    expr = expr[!is.na(rownames(expr)) & rownames(expr) != "",]

    genes = expr[genes,,drop=FALSE]
    expr = expr[rownames(expr) != "STAT1",]
    gsva = GSVA::gsva(expr, sets, parallel.sz=1)
    scores = rbind(genes, gsva)
}

# response signatures: Stat1 ChIP (enrichr), Ifn response
genes = c("STAT1", "PIAS1", "IFNG", "TP53")
ee = io$load("../expr_diff/eset_Mad2PB.RData")$vs
stat1 = ee["Stat1",]
rest = ee[-which(rownames(ee) == "Stat1"),]
cc = data.frame(
    gene = toupper(rownames(rest)),
    cor = narray::map(rest, along=2, function(x) cor(x, stat1))
) %>% arrange(-cor)

encc = io$load("../data/genesets/human/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.RData")
chea = io$load("../data/genesets/human/ChEA_2016.RData")
go = io$load("../data/genesets/human/GO_Biological_Process_2018.RData")
hm = io$load("../data/genesets/human/CH.HALLMARK.RData")
mile = io$load("../data/genesets/human/MILE_regulons.RData")
sets = c(go["regulation of complement activation (GO:0030449)"],
         chea["STAT1_20625510_ChIP-Seq_HELA_Human"],
         chea["STAT1_17558387_ChIP-Seq_HELA_Human"],
         encc["STAT3_CHEA"],
         encc["STAT3_ENCODE"],
         hm["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
         STAT1_complement = list(intersect(
            go[["regulation of complement activation (GO:0030449)"]],
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
         STAT1_mhc = list(intersect(
            go[["antigen receptor-mediated signaling pathway (GO:0050851)"]],
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
         STAT1_apop = list(intersect(
            go[["apoptotic process (GO:0006915)"]],
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
         STAT1_cor = list(intersect(
            head(cc$gene, 1000),
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]]))
)
sets$STAT1_IFNg = intersect(sets$HALLMARK_INTERFERON_GAMMA_RESPONSE,
                            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])
sets$noSTAT1_IFNg = setdiff(sets$HALLMARK_INTERFERON_GAMMA_RESPONSE,
                            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])
sets$HALLMARK_INTERFERON_GAMMA_RESPONSE = NULL
# + MHC, IFN(g), apop(?) -> do GO enrichment, then choose (or: clust?)

cohorts = c("BRCA", "LUAD", "COAD", "PRAD", "SKCM")
gsva = lapply(cohorts, gsva_tcga, sets=sets, genes=genes)
scores = gsva %>%
    narray::stack(along=2) %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample") %>%
    mutate(cohort = tcga$barcode2study(sample)) %>%
    filter(substr(sample, 14, 16) == "01A") %>%
    mutate(sample = substr(sample, 1,12))
#    tcga$filter(cancer==TRUE, primary==TRUE)

immune = readxl::read_xlsx("1-s2.0-S1074761318301213-mmc2.xlsx") %>%
    transmute(sample = `TCGA Participant Barcode`,
              aneup = as.numeric(`Aneuploidy Score`),
              stroma = as.numeric(`Stromal Fraction`),
              immune = as.numeric(Lymphocytes),
              NK_activated = as.numeric(`NK Cells Activated`),
              NK_total = NK_activated + as.numeric(`NK Cells Resting`))

data = inner_join(scores, immune, by="sample") %>% na.omit()
mat = t(data.matrix(data[,! colnames(data) %in% c("sample", "cohort")]))
colnames(mat) = data$sample
groups = narray::mask(data$cohort) + 0
tmat = narray::split(mat, along=2, subsets=data$cohort)
mat2 = rbind(mat, t(groups))

pdf("gsva.pdf", 20, 15)
util$plot_cor_matrix(t(mat2), text_color=NULL)
util$plot_cor_matrix(t(mat), text_color=NULL)

util$pcor(t(mat)) %>% util$plot_pcor_net(node_size=4, edge_size=2.5)
util$plot_bootstrapped_pcor(t(mat), node_size=4)
util$pcor(t(mat2)) %>% util$plot_pcor_net(node_size=4, edge_size=2.5)
util$plot_bootstrapped_pcor(t(mat2), node_size=4)

for (i in seq_along(tmat)) {
    name = names(tmat)[i]
    util$plot_cor_matrix(t(tmat[[i]]), title=name, text_color=NULL)
    util$plot_bootstrapped_pcor(t(tmat[[i]]), fdr=0.3, node_size=4, title=name)
}
dev.off()
