# test for insertion differences with aneuploidy
# plot the result as a volcano (dna, rna separate)/heatmap (both dna+rna)
library(dplyr)
io = import('io')
sys = import('sys')

#' Plot read counts for all samples/exons with anntotation insertion sites
plot_exon_expr = function() {
}

args = sys$cmd$parse(
    opt('i', 'insertions', 'all DNA insertions table', 'cis_sanger.tsv'),
    opt('d', 'dna', 'CIS in DNA', 'cis_sanger_results.tsv'),
    opt('r', 'rna', 'CTG in RNA', '../data/rnaseq_imfusion/merged_ctgs.txt'),
    opt('e', 'exons', 'exon expression table', '../data/rnaseq_imfusion/exon_counts.txt'),
    opt('p', 'plotfile', 'pdf to plot to', 'plot_aneup.pdf'))

ins = io$read_table(args$insertions)
expr = io$read_table(args$exons)
cis = io$read_table(args$dna)
ctg = io$read_table(args$rna)
