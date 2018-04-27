library(dplyr)
b = import('base')
io = import('io')

# load and merge aneuploidy scores
aneup = io$read_table("../ploidy_compare/compare_ploidy.tsv", header=TRUE) %>%
    select(-coverage, -tissue) %>%
    filter(!duplicated(data.frame(type, sample))) %>%
    tidyr::spread("type", "aneuploidy") %>%
    group_by(sample) %>%
    summarize(aneup = ifelse(sample == "",
        `WGS (single-cell)`,
        `WGS (merged)` %or% `WGS (30-cell)`)) %>%
    arrange(aneup) %>%
    filter(!sample %in% c("S", "T")) %>%
    mutate(sample = paste0(sub("[ST]", "", sample), tolower(sub("[0-9]+", "", sample))))

cis = io$load("../data/cis/cis_per_tumor.RData") %>%
    mutate(sample = ifelse(sample == "410bm", "410t", sample),
           sample = ifelse(sample == "415bm", "415s", sample),
           sample = ifelse(sample == "452s", "452t", sample)) %>%
    filter(sample != "477s") %>%
    inner_join(aneup)

save(cis, file="dset.RData")
