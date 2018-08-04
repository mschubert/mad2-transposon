library(DESeq2)
library(ggplot2)
library(dplyr)
io = import('io')
sys = import('sys')
util = import('./aneup_de')

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'gene expression RData', 'aneup_de.RData'),
    opt('i', 'ins', 'gene name of insert', 'Erg'),
    opt('o', 'outfile', 'results RData', 'aneup_de_Erg.RData'),
    opt('p', 'plotfile', 'pdf', 'aneup_de_Erg.pdf'))

dset = io$load(args$diff_expr)
cis = dset$cis$samples %>%
    filter(external_gene_name == args$ins) %>%
    select(sample, ins=external_gene_name)
eset = dset$eset
idx = colData(eset) %>%
    as.data.frame() %>%
    left_join(cis) %>%
    mutate(ins = ifelse(is.na(ins), 0, 1),
           noins = (!ins) + 0)
eset@colData = DataFrame(idx)

design(eset) = ~ tissue + type + aneup * ins
res = DESeq2::DESeq(eset)
DESeq2::resultsNames(res)

pdf(args$plotfile)
print(util$plot_pcs(idx, eset$pca, 1, 2)) #todo: annotate ins
for (name in c("aneup", "ins", "aneup.ins"))
    print(util$plot_volcano(res, name))
dev.off()

save(res, file=args$outfile)
