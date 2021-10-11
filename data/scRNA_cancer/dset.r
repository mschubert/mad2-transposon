library(dplyr)

bpath = "/home/p.batistatan/main/NKI-CIN-internship-main-backup"

dsets = c(
    "BakhoumEtAl2018",
    "ChungEtAl2017",
    "DarmanisEtAl2017",
    "GiustacchiniEtAl2017",
    "KaraayvazEtAl2018",
    "LeeEtAl2020",
    "NelsonEtAl2020",
    "PuramEtAl2017",
    "SunEtAl2020",
    "TijhuisEtAl202X",
    "TiroshEtAl2016"
)

res = file.path(bpath, dsets, "output_DE/df.deseq.results.pseudobulkWald.csv") %>%
    lapply(read.table, sep=";", dec=",", header=TRUE) %>%
    lapply(as_tibble) %>%
    setNames(dsets)

saveRDS(res, file="dset.rds")
