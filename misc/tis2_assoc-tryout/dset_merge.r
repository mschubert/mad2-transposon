library(dplyr)
b = import('base')
io = import('io')

# source sample IDs are both 123S and T567
fix_ids = function(x) {
    letter = tolower(sub("[0-9]+", "", x))
    nums = sub("[ST]", "", x)
    paste0(nums, letter)
}

# load and merge aneuploidy scores
aneup = io$read_table("../ploidy_compare/compare_ploidy.tsv", header=TRUE) %>%
    select(-coverage, -tissue) %>%
    filter(!duplicated(data.frame(type, sample)),
           !sample %in% c("S", "T")) %>%
    mutate(sample = fix_ids(sample)) %>%
    tidyr::spread("type", "aneuploidy") %>%
    group_by(sample) %>%
    summarize(aneup = ifelse(sample == "",
        `WGS (single-cell)`,
        `WGS (merged)` %or% `WGS (30-cell)`)) %>%
    arrange(-aneup)

cis = io$load("../data/cis/cis_per_tumor.rds") %>%
    mutate(sample = ifelse(sample == "410bm", "410t", sample),
           sample = ifelse(sample == "415bm", "415s", sample)) %>%
    filter(sample != "477s") %>%
    inner_join(aneup)

save(cis, file="dset.RData")
