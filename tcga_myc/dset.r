library(dplyr)
library(ggplot2)
library(patchwork)
io = import('io')
sys = import('sys')
tcga = import('data/tcga') # todo: include this in dset

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', 'genesets.yaml'),
    opt('t', 'tcga', 'rds', 'tcga.rds'),
    opt('o', 'outfile', 'rds', 'dset.rds'),
    arg('setfiles', '', arity='*',
        list.files("gsva", "\\.rds$", recursive=TRUE, full.names=TRUE))
)

cfg = yaml::read_yaml(args$config)
sets = lapply(args$setfiles, readRDS) %>%
    setNames(tools::file_path_sans_ext(basename(args$setfiles)))
#colnames(sets$SPEED2) = paste0("speed2.", colnames(sets$SPEED2))

# process gene set data
#prog = io$load("../../speed2/scores/tcga/pathways_mapped/speed_matrix.RData")
#prog = prog[tcga$barcode2study(rownames(prog)) == "LIHC",]
#colnames(prog) = paste0("prog.", colnames(prog))

immdf = readxl::read_xlsx("1-s2.0-S1074761318301213-mmc2.xlsx", 1) %>%
    dplyr::rename(Sample = `TCGA Participant Barcode`,
                  Study = `TCGA Study`) %>%
    filter(Study == "BRCA")
immune = data.matrix(immdf[-c(1:4)])
rownames(immune) = paste0(immdf$Sample, "-01A")

smat = narray::stack(c(sets, list(immune)), along=2)
subs = unlist(cfg$dset[-1], use.names=FALSE)
smat = smat[,subs]

# process TCGA data
dset = readRDS(args$tcga)
meta = dset$meta
expr = t(dset$expr)
cna = ((t(dset$cna) - 2) / meta$purity) + 2

mut = function(g) meta$Sample %in% (filter(dset$mut, Hugo_Symbol %in% g) %>% pull(Tumor_Sample_Barcode))
amp = function(g) cna[,g] > 2.5
del = function(g) cna[,g] < 1.5

narray::intersect(meta$Sample, expr, cna, smat, along=1)
dmat = cbind(p53_mut = mut("TP53"), # | del("TP53"),
             p53_copy = cna[,"TP53"],
             pi3k_mut = mut(c("PIK3CA", "PTEN")) | del("PTEN"),
             pten_copy = cna[,"PTEN"],
             mapk_mut = mut(c("KRAS", "HRAS", "NRAS")),
             myc_amp = amp("MYC"),
             myc_copy = cna[,"MYC"], smat)
#             expr[,c("IL1RN","MYC")], smat)

saveRDS(list(dmat=dmat, smat=smat, meta=meta), file=args$outfile)
