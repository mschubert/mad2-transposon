library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'all samples .RData', '../data/cis/cis_per_tumor.RData'),
    opt('r', 'reads', 'min reads to consider ins', '20'),
    opt('o', 'outfile', 'filtered samples & positions', 'analysis_set.RData'))

#TODO: read this from sample sheet (?)
exclude = c(
    "401s", # using t
    "408s", # ctl
    "410s", # only one transposon; use uLN instead
    "414s", # ctl
    "418s", # using t
    "422t", # using s
    "426s", # ctl
    "430s", # ctl
    "436s", # ctl
    "438s", # ctl
    "439s", # ctl
    "441s", # ctl
    "445s", # ctl
    "446t", # using s
    "452s", # using t
    "460s", # ctl
    "475s", # ctl
    "477s", # using t
    "478s", # control
    "480s", # control
    "610s",
    "612s", # using t
    "613s", # using t
    NULL
)

replace = c(
    "1jspbpb" = "145s",
    "2jspbpb" = "155s",
    "3jspbpb" = "157s",
    "4jspbpb" = "180s",
    "8jspbpb" = "184s",
    "23jspbpb" = "212s",
    "410uln" = "410s",
    NULL
)

# s&t different tumors: 184, 443 (RNA also both)
ins = io$load(args$infile) %>%
    dplyr::select(sample, chr, position, reads) %>% # ignore strand
    distinct() %>%
    filter(!sample %in% exclude) %>%
    mutate(sample = ifelse(sample %in% names(replace), replace[sample], sample)) %>%
    filter(grepl("[0-9]{3}[st]", sample),
           reads >= as.integer(args$reads)) %>%
    arrange(sample, chr, position)

save(ins, file=args$outfile)
