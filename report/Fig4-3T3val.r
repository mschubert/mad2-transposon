library(dplyr)
library(ggplot2)
sys = import('sys')

read_pzfx = function(fname) {
    read_sheet = function(sheet) {
        pzfx::read_pzfx(fname, sheet) %>%
            as_tibble(.name_repair="universal") %>%
            tidyr::gather("sample", "weight", -Var.1) %>%
            filter(!is.na(weight)) %>%
            mutate(genotype = sub(".*(myc|p53).*", "\\1", sample),
                   CIN = sub(".*(dnMCAK|KIF2C).*", "\\1", sample),
                   STAT1 = ifelse(grepl("STAT1\\.KO", sample), "ko", "wt"),
                   mouse = paste("myc", CIN, STAT1, sub(".*([0-9]+)$", "\\1", sample))) %>%
            select(genotype, CIN, STAT1, mouse, day=Var.1, weight)
    }
    lapply(pzfx::pzfx_tables(fname), read_sheet) %>%
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

#    pzfx::read_pzfx("external/3t3/Balbc both replicates immune landscape.pzfx")
})
