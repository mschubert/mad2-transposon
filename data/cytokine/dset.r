library(dplyr)
library(ggplot2)
plt = import('plot')
sys = import('sys')
theme_set(cowplot::theme_cowplot())

args = sys$cmd$parse(
    opt('i', 'infile', 'xlsx with results', 'PR-080-TES_RESULTS REPORT.xlsx'),
    opt('o', 'outfile', 'rds', 'dset.rds'),
    opt('p', 'plotfile', 'pdf', 'dset.pdf'))

### phospho measurements ###
pdata = readxl::read_xlsx(args$infile, "RESULTS-PLATE ID PL1004", range="B17:S91")
pmat = data.matrix(pdata[-1])
rownames(pmat) = pdata$SAMPLE
annot1 = data.frame(SAMPLE=rownames(pmat)) %>%
    mutate(time = as.integer(sub(".*T=([0-9]+)$", "\\1", SAMPLE)),
           cells = case_when(grepl("BT549", SAMPLE) ~ "BT549",
                             grepl("RPE", SAMPLE) ~ "RPE-1",
                             TRUE ~ ""),
           treatment = case_when(grepl("DMSO", SAMPLE) ~ "DMSO",
                                 grepl("rev", SAMPLE) ~ "rev",
                                 grepl("IFN", SAMPLE) ~ "IFNg",
                                 TRUE ~ ""))
phospho = pdata %>%
    tidyr::gather("protein", "value", -SAMPLE) %>%
    left_join(annot1, by="SAMPLE")

ridge1 = phospho %>%
    filter(cells != "" & !is.na(time)) %>%
    mutate(value = pmin(value, 10000)) %>%
    ggplot(aes(x=value, y=protein, fill=treatment)) +
        ggridges::geom_density_ridges(alpha=0.4) +
        facet_grid(cells ~ time) +
        scale_fill_manual(values=c("lightblue", "cyan", "red")) +
        ggtitle("Phospho measurements")

pca1 = prcomp(pmat, scale.=TRUE)
p11 = plt$pca(pca1, aes(x=PC1, y=PC2), annot=annot1, biplot=TRUE) +
    geom_point(aes(shape=cells, size=time, color=treatment)) +
    geom_text_repel(aes(label=SAMPLE), size=2)
p12 = plt$pca(pca1, aes(x=PC3, y=PC4), annot=annot1, biplot=TRUE) +
    geom_point(aes(shape=cells, size=time, color=treatment)) +
    geom_text_repel(aes(label=SAMPLE), size=2)

### cytokine measurements ###
cdata = readxl::read_xlsx(args$infile, "RESULTS-PLATE ID PL1005", range="A29:U119") %>%
    filter(!grepl("Standard", SAMPLE))
cmat = data.matrix(cdata[-1])
rownames(cmat) = cdata$SAMPLE
annot2 = data.frame(SAMPLE=rownames(cmat)) %>%
    mutate(time = as.integer(sub(".*T=([0-9]+)$", "\\1", SAMPLE)),
           cells = case_when(grepl("BT549", SAMPLE) ~ "BT549",
                             grepl("RPE", SAMPLE) ~ "RPE-1",
                             TRUE ~ ""),
           treatment = case_when(grepl("DMSO", SAMPLE) ~ "DMSO",
                                 grepl("rev", SAMPLE) ~ "rev",
                                 grepl("IFN", SAMPLE) ~ "IFNg",
                                 TRUE ~ ""))
cytokine = cdata %>%
    tidyr::gather("protein", "value", -SAMPLE) %>%
    group_by(protein) %>%
        mutate(value = value - mean(value[SAMPLE == "Background"])) %>%
    ungroup() %>%
    filter(SAMPLE != "Background") %>%
    left_join(annot2, by="SAMPLE")

ridge2 = cytokine %>%
    filter(cells != "" & !is.na(time)) %>%
    group_by(protein) %>%
        mutate(z_value = scale(value)) %>%
    ungroup() %>%
    ggplot(aes(x=z_value, y=protein, fill=treatment)) +
        ggridges::geom_density_ridges(alpha=0.4) +
        facet_grid(cells ~ time) +
        scale_fill_manual(values=c("lightblue", "cyan", "red")) +
        ggtitle("Cytokine measurements")

pca2 = prcomp(cmat, scale.=TRUE)
p21 = plt$pca(pca2, aes(x=PC1, y=PC2), annot=annot2, biplot=TRUE) +
    geom_point(aes(shape=cells, size=time, color=treatment)) +
    geom_text_repel(aes(label=SAMPLE), size=2)
p22 = plt$pca(pca2, aes(x=PC3, y=PC4), annot=annot2, biplot=TRUE) +
    geom_point(aes(shape=cells, size=time, color=treatment)) +
    geom_text_repel(aes(label=SAMPLE), size=2)

pdf(args$plotfile, 16, 12)
print(ridge1)
print(p11)
print(p12)
print(ridge2)
print(p21)
print(p22)
dev.off()

saveRDS(cytokine=cytokine, phospho=phopho, file=args$outfile)
