library(dplyr)
library(ggplot2)
library(patchwork)
io = import('io')
sys = import('sys')
seq = import('seq')
tcga = import('data/tcga') # todo: include this in dset

read_sets = function(setname) {
    sets = readRDS(grep(setname, args$setfiles, value=TRUE))
    t(sets[cfg$dset[[setname]],,drop=FALSE])
}

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', 'genesets.yaml'),
    opt('t', 'tcga', 'rds', 'tcga.rds'),
    opt('o', 'outfile', 'rds', 'dset.rds'),
    arg('setfiles', '', arity='*',
        list.files("../data/gsva/tcga-brca", "\\.rds$", recursive=TRUE, full.names=TRUE))
)

cfg = yaml::read_yaml(args$config)
sets = setdiff(names(cfg$dset), "genes") %>%
    sapply(read_sets, simplify=FALSE)

immune = tcga$cell_frac("BRCA", "epic")

smat = narray::stack(c(sets, list(immune)), along=2)

# process TCGA data
dset = readRDS(args$tcga)
meta = dset$meta %>%
    inner_join(tcga$aneuploidy("BRCA")) %>%
    mutate(aneup_log2seg = aneup_log2seg / purity)
expr = t(dset$expr)
cna = ((t(dset$cna) - 2) / meta$purity) + 2

mut = function(g) meta$Sample %in% (filter(dset$mut, Hugo_Symbol %in% g) %>% pull(Tumor_Sample_Barcode))
amp = function(g) cna[,g] > 2.5
del = function(g) cna[,g] < 1.5

narray::intersect(meta$Sample, expr, cna, smat, along=1)
dmat = cbind(#ifn_loss = del("IFNA1"),
             p53_mut = mut("TP53"), # | del("TP53"),
             p53_copy = cna[,"TP53"],
             pi3k_mut = mut(c("PIK3CA", "PTEN")) | del("PTEN"),
             pten_copy = cna[,"PTEN"],
             mapk_mut = mut(c("KRAS", "HRAS", "NRAS")),
             myc_amp = amp("MYC"),
             myc_copy = cna[,"MYC"],
             expr[,cfg$dset$genes$expr], smat)

saveRDS(list(dmat=dmat, smat=smat, meta=meta), file=args$outfile)
