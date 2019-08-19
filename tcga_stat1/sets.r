library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'sets.rds'))

# response signatures: Stat1 ChIP (enrichr), Ifn response
#ee = io$load("../expr_diff/eset_Mad2PB.RData")$vs
#stat1 = ee["Stat1",]
#rest = ee[-which(rownames(ee) == "Stat1"),]
#cc = data.frame(
#    gene = toupper(rownames(rest)),
#    cor = narray::map(rest, along=2, function(x) cor(x, stat1))
#) %>% arrange(-cor)

encc = io$load("../data/genesets/human/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.RData")
chea = io$load("../data/genesets/human/ChEA_2016.RData")
go = io$load("../data/genesets/human/GO_Biological_Process_2018.RData")
hm = io$load("../data/genesets/human/CH.HALLMARK.RData")
msdb = io$load("../data/genesets/human/MSigDB_Oncogenic_Signatures.RData")
reac = io$load("../data/genesets/human/Reactome_2016.RData")
#mile = io$load("../data/genesets/human/MILE_regulons.RData")
sets = c(go["regulation of complement activation (GO:0030449)"],
         chea["STAT1_20625510_ChIP-Seq_HELA_Human"],
         chea["STAT1_17558387_ChIP-Seq_HELA_Human"],
         encc["STAT3_CHEA"],
         encc["STAT3_ENCODE"],
         hm["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
         hm["HALLMARK_MYC_TARGETS_V2"],
         hm["HALLMARK_OXIDATIVE_PHOSPHORYLATION"],
         msdb["TBK1.DF_DN"],
         reac["Mitochondrial translation_Homo sapiens_R-HSA-5368287"],
         reac["Interferon alpha/beta signaling_Homo sapiens_R−HSA−909733"],
         go["mitochondrial respiratory chain complex assembly (GO:0033108)"],
#         go["type I interferon signaling pathway (GO:0060337)"], # kills IFNg>stat1>IFNg response
         STAT1_complement = list(intersect(
            go[["regulation of complement activation (GO:0030449)"]],
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
         STAT1_mhc = list(intersect(
            go[["antigen receptor-mediated signaling pathway (GO:0050851)"]],
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]])),
         STAT1_apop = list(intersect(
            go[["apoptotic process (GO:0006915)"]],
            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]]))
#         STAT1_cor = list(intersect(
#            head(cc$gene, 1000),
#            chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]]))
)
#sets$STAT1_IFNg = intersect(sets$HALLMARK_INTERFERON_GAMMA_RESPONSE,
#                            intersect(chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]],
#                                      chea["STAT1_20625510_ChIP-Seq_HELA_Human"]))
#sets$noSTAT1_IFNg = setdiff(sets$HALLMARK_INTERFERON_GAMMA_RESPONSE,
#                            c(chea[["STAT1_17558387_ChIP-Seq_HELA_Human"]],
#                              chea["STAT1_20625510_ChIP-Seq_HELA_Human"]))
#
#sets$HALLMARK_INTERFERON_GAMMA_RESPONSE = NULL
## + MHC, IFN(g), apop(?) -> do GO enrichment, then choose (or: clust?)

saveRDS(sets, file=args$outfile)
