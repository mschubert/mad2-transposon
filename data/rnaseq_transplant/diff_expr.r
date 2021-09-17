library(dplyr)
library(ggplot2)
sys = import('sys')

args = sys$cmd$parse(
    opt('e', 'eset', 'rds', 'samples.rds'),
    opt('c', 'comps', 'yaml', 'diff_expr.yaml'),
    opt('o', 'outfile', 'rds', 'diff_expr.rds'),
    opt('p', 'plotfile', 'pdf', 'diff_expr.pdf')
)
