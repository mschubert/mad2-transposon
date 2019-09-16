library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('o', 'outfile', 'rds', 'sets.rds'))

encc = io$load("../data/genesets/human/ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.RData")
chea = io$load("../data/genesets/human/ChEA_2016.RData")
go = io$load("../data/genesets/human/GO_Biological_Process_2018.RData")
hm = io$load("../data/genesets/human/CH.HALLMARK.RData")
msdb = io$load("../data/genesets/human/MSigDB_Oncogenic_Signatures.RData")
reac = io$load("../data/genesets/human/Reactome_2016.RData")
#mile = io$load("../data/genesets/human/MILE_regulons.RData")
sets1 = c(go["regulation of complement activation (GO:0030449)"],
          chea["STAT1_20625510_ChIP-Seq_HELA_Human"],
          chea["STAT1_17558387_ChIP-Seq_HELA_Human"],
          encc["STAT3_CHEA"],
          encc["STAT3_ENCODE"]
)
sets2 = c(hm["HALLMARK_INTERFERON_GAMMA_RESPONSE"],
         reac["Interferon alpha/beta signaling_Homo sapiens_R−HSA−909733"]
)
genes = intersect(unlist(sets1, use.names=FALSE),
                  unlist(sets2, use.names=FALSE))

saveRDS(list(sets=list(), genes=genes), file=args$outfile)
