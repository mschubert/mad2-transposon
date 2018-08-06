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
           noins = (!ins) + 0,
           aneup_ins = 2*(ins-0.5) * aneup)
eset@colData = DataFrame(idx)

design(eset) = ~ tissue + type + ins * aneup
res1 = DESeq2::estimateDispersions(eset) %>%
    DESeq2::nbinomWaldTest(maxit=1000)

design(eset) = ~ tissue + type + ins + aneup_ins
res2 = DESeq2::estimateDispersions(eset) %>%
    DESeq2::nbinomWaldTest(maxit=1000)

pdf(args$plotfile)
print(util$plot_pcs(idx, dset$pca, 1, 2, hl=cis$sample))
for (name in c("aneup", "ins", "ins.aneup"))
    print(util$plot_volcano(res1, name))
print(util$plot_volcano(res2, "aneup_ins"))
dev.off()

save(res1, res2, file=args$outfile)
