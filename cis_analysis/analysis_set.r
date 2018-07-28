library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'all samples .RData', '../data/cis/cis_per_tumor.RData'),
    opt('r', 'reads', 'min reads to consider ins', '20'),
    opt('o', 'outfile', 'filtered samples & positions', 'analysis_set.RData'))

#TODO: read this from sample sheet
exclude = c(
    "401s", # using t
    "408s", # ctl
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
    "621s", # ???
    "628s", # ???
    "629s", # ???
    NULL
)

# s&t different tumors: 184, 443 (RNA also both)
# missing: 145s, 155s, 157s, 180s, 184s, 410s

ins = io$load(args$infile) %>%
    dplyr::select(sample, chr, position, reads) %>% # ignore strand
    distinct() %>%
    filter(!sample %in% exclude,
           grepl("[0-9]{3}[st]", sample),
           reads >= as.integer(args$reads)) %>%
    arrange(sample, chr, position)

save(ins, file=args$outfile)
