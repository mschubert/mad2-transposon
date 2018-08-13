library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'xls read', 'Final table mice.xlsx'),
    opt('o', 'outfile', 'tsv', 'meta.tsv'))

tab = readxl::read_excel(args$infile, na="n.a.") %>%
    filter(!is.na(Mouse)) %>%
    transmute(mouse = Mouse,
              genotype = Genotype,
              sex = Gender,
              months_injection = `Time after poly (I:C) (m)`,
              months_death = `Age at death (m)`,
              tissue = tolower(`Primary affected organ`),
              hist_nr = `Hist nr.`,
              sample = paste0(hist_nr, substr(tissue, 0, 1)),
              wbc_per_nl = as.numeric(ifelse(`WBC (·109/L)` == ">100", 100, `WBC (·109/L)`)),
              spleen_g = `Spleen (mg)` / 1000,
              thymus_g = `Thymus (mg)` / 1000,
              type = factor(sub("^([^ ]+).*", "\\1", Pathology), levels=c("Myeloid", "T-cell", "Other")))

io$write_table(tab, file=args$outfile)
