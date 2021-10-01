library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'sets.rds')
)

genes = c("STAT1", "PIAS1", "IFNG", "IL1B", "IFITM1",
          "CGAS", "TBK1", "IRF2", "IRF3", "TP53")
# don't have: "IFNA", "IFNB", "OAS1", "ISG54"

encc = readRDS("../data/genesets/human/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.rds")
chea = readRDS("../data/genesets/human/ChEA_2016.rds")
go = readRDS("../data/genesets/human/GO_Biological_Process_2021.rds")
hm = readRDS("../data/genesets/human/MSigDB_Hallmark_2020.rds")
msdb = readRDS("../data/genesets/human/MSigDB_Oncogenic_Signatures.rds")
reac = readRDS("../data/genesets/human/Reactome_2016.rds")
sets = c(go["regulation of complement activation (GO:0030449)"],
         chea["STAT1 20625510 ChIP-Seq HELA Human"],
         chea["STAT1 17558387 ChIP-Seq HELA Human"],
         encc["STAT3 CHEA"],
         encc["STAT3 ENCODE"],
         hm["Interferon Gamma Response"],
         hm["Myc Targets V2"],
         hm["Oxidative Phosphorylation"],
         msdb["TBK1.DF DN"],
         reac["Mitochondrial translation Homo sapiens R-HSA-5368287"],
         reac["Interferon Signaling Homo sapiens R-HSA-913531"],
         go["mitochondrial respiratory chain complex assembly (GO:0033108)"],
#         go["type I interferon signaling pathway (GO:0060337)"], # kills IFNg>stat1>IFNg response
         STAT1_complement = list(intersect(
            go[["regulation of complement activation (GO:0030449)"]],
            chea[["STAT1 17558387 ChIP-Seq HELA Human"]])),
         STAT1_mhc = list(intersect(
            go[["antigen receptor-mediated signaling pathway (GO:0050851)"]],
            chea[["STAT1 17558387 ChIP-Seq HELA Human"]])),
         STAT1_apop = list(intersect(
            go[["apoptotic process (GO:0006915)"]],
            chea[["STAT1 17558387 ChIP-Seq HELA Human"]])),
         go["regulation of complement activation (GO:0030449)"],
         go["antigen receptor-mediated signaling pathway (GO:0050851)"],
         go["apoptotic process (GO:0006915)"]
)

genes = c(genes, sets$STAT1_mhc)

saveRDS(list(sets=sets, genes=genes), file=args$outfile)
