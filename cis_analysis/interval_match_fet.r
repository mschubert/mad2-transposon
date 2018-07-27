io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('d', 'data', 'cis RData', '../data/cis/cis_per_tumor.RData'),
    opt('o', 'outfile', 'insertion statistics RData', 'interval_match_fet.RData'))
