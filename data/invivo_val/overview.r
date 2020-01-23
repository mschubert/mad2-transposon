library(dplyr)
library(ggplot2)

clusters = readxl::read_xlsx("metaclusters_christy.xlsx")
cluster_annot = readr::read_tsv("clusters.txt")
mouse_annot = readxl::read_xlsx("tumor list.xlsx") %>%
    mutate(mouse = sub("\\ .*$", "", Tumor))
totals = readxl::read_xlsx("immune_cells_total.xlsx")

dset = tidyr::gather(clusters, "Tumor", "value", -MetaclusterID) %>%
    left_join(cluster_annot) %>%
    left_join(mouse_annot) %>%
    left_join(totals) %>%
    select(-Tumor, -MetaclusterID) %>%
    mutate(cGAS = paste("cGAS", cGAS),
           STAT1 = paste("STAT1", STAT1)) %>%
    rename(dnMCAK_over_WT = `WT/dnMCAK`)
dset = dset %>%
    group_by(mouse, cGAS, STAT1, dnMCAK_over_WT) %>%
    summarize(value = sum(value) / 2,
              `% immune cells` = unique(`% immune cells`),
              `immune cell count` = unique(`immune cell count`)) %>%
    mutate(cells = "total_immune_half") %>%
    bind_rows(dset)

dnMCAK_over_WT = dset %>%
    group_by(cells, cGAS, STAT1, mouse) %>%
    summarize(value = value[dnMCAK_over_WT == "dnMCAK"] / value[dnMCAK_over_WT == "WT"]) %>%
    ungroup()

abundance = dset %>% mutate(value=value * `% immune cells`)

ab_dnMCAK_over_WT = abundance %>%
    group_by(cells, cGAS, STAT1, mouse) %>%
    summarize(value = value[dnMCAK_over_WT == "dnMCAK"] / value[dnMCAK_over_WT == "WT"]) %>%
    ungroup()

p1 = ggplot(abundance, aes(x=cells, y=value, color=dnMCAK_over_WT)) +
    ggbeeswarm::geom_quasirandom() +
    facet_wrap( ~ cGAS+STAT1) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=0.5)) +
    ggtitle("effect of genotypes on immune cell abundance")

p2 = ggplot(ab_dnMCAK_over_WT, aes(x=cells, y=log2(value), color=STAT1)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point() +
    facet_wrap( ~ cGAS) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=0.5)) +
    ggtitle("effect of dnMCAK over WT on immune celabundance")

p3 = ggplot(dset, aes(x=cells, y=value, color=dnMCAK_over_WT)) +
    ggbeeswarm::geom_quasirandom() +
    facet_wrap( ~ cGAS+STAT1) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=0.5)) +
    ggtitle("effect of genotypes on immune cell composition")

p4 = ggplot(dnMCAK_over_WT, aes(x=cells, y=log2(value), color=STAT1)) +
    geom_boxplot(outlier.shape=NA) +
    geom_point() +
    facet_wrap( ~ cGAS) +
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=0.5)) +
    ggtitle("effect of dnMCAK over WT on immune cell composition")

pdf("overview.pdf", 12, 10)
print(p1)
print(p2)
print(p3)
print(p4)
dev.off()
