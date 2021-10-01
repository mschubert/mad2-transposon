library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'sets.rds')
)

encc = readRDS("../data/genesets/human/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.rds")
chea = readRDS("../data/genesets/human/ChEA_2016.rds")
go = readRDS("../data/genesets/human/GO_Biological_Process_2021.rds")
hm = readRDS("../data/genesets/human/MSigDB_Hallmark_2020.rds")
msdb = readRDS("../data/genesets/human/MSigDB_Oncogenic_Signatures.rds")
reac = readRDS("../data/genesets/human/Reactome_2016.rds")
#mile = readRDS("../data/genesets/human/MILE_regulons.rds")
sets1 = c(go["regulation of complement activation (GO:0030449)"],
          chea["STAT1 20625510 ChIP-Seq HELA Human"],
          chea["STAT1 17558387 ChIP-Seq HELA Human"],
          encc["STAT3 CHEA"],
          encc["STAT3 ENCODE"]
)
sets2 = c(hm["Interferon Gamma Response"],
         reac["Interferon alpha/beta signaling Homo sapiens R-HSA-909733"]
)
genes = intersect(unlist(sets1, use.names=FALSE),
                  unlist(sets2, use.names=FALSE))

saveRDS(list(sets=list(), genes=genes), file=args$outfile)
