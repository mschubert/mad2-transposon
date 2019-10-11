library(ggplot2)
theme_set(cowplot::theme_cowplot())
gdsc = import('data/gdsc')

keep = gdsc$tissues(c("LAML", "DLBC", "LCML", "ALL"))
expr = t(gdsc$basal_expression())[,c("ETS1", "ETS2", "ERG")]
aneup = gdsc$aneuploidy()
narray::intersect(keep, expr, aneup$cosmic, along=1)

df = data.frame(name=aneup$cell_line, tissue=keep, aneuploidy=aneup$aneuploidy, expr)
ggplot(df, aes(x=ETS1, y=ERG)) +
    geom_point(aes(shape=tissue, color=aneuploidy), size=5) +
    geom_text_repel(aes(label=name), size=2)

ggplot(df %>% filter(keep == "ALL"), aes(x=ETS1, y=ERG)) +
    geom_point(aes(shape=tissue, color=aneuploidy), size=5) +
    geom_text_repel(aes(label=name), size=2)

dev.off()
