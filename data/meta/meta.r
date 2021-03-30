library(dplyr)
#library(googlesheets4) # why do we need auth?
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'xls read', 'TPS_screen_MasterList.xlsx'),
    opt('a', 'aneup', '../ploidy_compare/analysis_set.rds', NA),
    opt('o', 'outfile', 'rds', 'meta.rds')
)

tab = readxl::read_excel(args$infile, na="n.a.") %>%
    transmute(mouse = `Mouse no`,
              genotype = Genotype,
              sex = Sex,
              months_injection = `Time after poly (I:C) (m)`,
              months_death = `Age at death (m)`,
              tissue = tolower(`Primary affected organ`),
              hist_nr = `Hist nr.`,
              sample = `Sample ID`,
#              wbc_per_nl = as.numeric(ifelse(`WBC (1e9/L)` == ">100", 100, `WBC (1e9/L)`)),
              spleen_g = `Spleen (mg)` / 1000,
              thymus_g = `Thymus (mg)` / 1000,
              type = factor(`Consensus type`, levels=c("Myeloid", "B-like", "T-cell")),
              subtype = factor(Subtype, levels=c("Ebf1", "Ets1", "Erg")),
              analysis_set = !is.na(`Analysis set`) & `Analysis set` == "y")

if (!is.na(args$aneup)) {
    aneup = readRDS(args$aneup)
    # add aneuploidy
}

saveRDS(meta, file=args$outfile)
