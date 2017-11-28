library(cowplot)
library(ggridges)
io = import('io')
idmap = import('process/idmap')

# reads per chromosome for T-ALLs
dset = io$load('../data/rnaseq/assemble.RData')
chrs = idmap$gene(rownames(dset$counts), to="chromosome_name", dset="mmusculus_gene_ensembl")
reads = dset$counts[!is.na(chrs),]
chrs = chrs[!is.na(chrs)]
reads = narray::map(reads, along=1, sum, subsets=chrs)[as.character(1:19),]

# reads per chromosome for euploid controls
ctl = readr::read_tsv("../../aneuploidy/data/rnaseq/T-ALL_read_count.txt")
ctl_reads = data.matrix(ctl[,c("eT_p0","eT_p2")])
chrs = idmap$gene(ctl$ensembl_gene_id, to="chromosome_name", dset="mmusculus_gene_ensembl")
ctl_reads = ctl_reads[!is.na(chrs),]
chrs = chrs[!is.na(chrs)]
ctl_reads = narray::map(ctl_reads, along=1, sum, subsets=chrs)[as.character(1:19),]

ploidy_multiples = reads / narray::crep(rowSums(ctl_reads), ncol(reads))
medians = narray::map(ploidy_multiples, along=1, median)
ploidy = 2 * ploidy_multiples / narray::rrep(medians, nrow(ploidy_multiples))

# compare to measured ploidy using scWGS
