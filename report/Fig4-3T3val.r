library(dplyr)
library(ggplot2)
sys = import('sys')

read_pzfx = function(fname, tables=pzfx::pzfx_tables(fname)) {
    read_sheet = function(sheet) {
        re = pzfx::read_pzfx(fname, sheet)
        colnames(re)[colnames(re) %in% c("Var.1", "ROWTITLE", "Day after transplantation")] = "day"
        re %>%
            as_tibble(.name_repair="universal") %>%
            tidyr::gather("sample", "weight", -day) %>%
            filter(!is.na(weight)) %>%
            mutate(genotype = sub(".*(myc|p53).*", "\\1", tolower(paste(sheet, sample))),
                   CIN = sub(".*(dnMCAK|Kif2c).*", "\\1", sub("KIF2C", "Kif2c", sample)),
                   STAT1 = ifelse(grepl("STAT1", sample), "ko", "wt"),
                   mouse = paste("myc", CIN, STAT1, sub(".*([0-9]+)$", "\\1", tolower(sample)))) %>%
            select(genotype, CIN, STAT1, mouse, day, weight)
    }
    lapply(tables, read_sheet) %>%
        bind_rows()
}

sys$run({
    args = sys$cmd$parse(
        opt('i', 'infile', 'prism pzfx', ''),
        opt('p', 'plotfile', 'pdf', 'Fig4-3T3val.pdf')
    )

    athy = read_pzfx("external/3t3/Growth curve athymic nude - Sept 2021_v2.pzfx")
    ggplot(athy, aes(x=day, y=weight, group=mouse)) +
        geom_line(aes(color=CIN)) +
        facet_grid(genotype ~ STAT1)

    balb = list(
        rep1 = read_pzfx("external/3t3/Fig5- Tumor plot Balbc rep 1.pzfx", c("Myc group", "p53 group")) %>%
            mutate(mouse = paste(mouse, "rep1")),
        rep2 = read_pzfx("external/3t3/Fig5- Tumor plot Balbc rep 2.pzfx", c("myc oex", "p53 KO")) %>%
            mutate(mouse = paste(mouse, "rep"))
    ) %>% bind_rows(.id="rep")
    ggplot(balb, aes(x=day, y=weight, group=mouse)) +
        geom_line(aes(color=CIN)) +
        facet_grid(genotype ~ STAT1)
})
