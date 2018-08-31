library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'cis RData', '../data/cis/cis_per_tumor.RData'),
    opt('o', 'outfile', 'tsv to save to', 'cis_sanger.tsv'))

hits = io$load(args$infile) %>%
    filter(hit_dist == 0) %>%
    group_by(sample) %>%
    mutate(i = seq_len(n())) %>%
    ungroup() %>%
    transmute(sample = sample,
              id = sprintf("%s.INS_%i", sample, i),
              seqname = chr,
              position = position,
              strand = ifelse(strand == "+", 1, -1),
              support = reads,
              support_spanning = 1,
              support_junction = "DNA",
              gene_id = ensembl_gene_id,
              gene_name = gene_name)

io$write_table(hits, file=args$outfile)
