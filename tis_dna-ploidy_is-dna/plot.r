library(dplyr)
library(patchwork)
b = import('base')
io = import('io')
plt = import('plot')
sys = import('sys')

args = sys$cmd$parse(
    opt('a', 'assocs', 'assocs .RData', 'betareg-1clonal.RData'),
    opt('s', 'subset', 'assoc object o use', 'hits'),
    opt('p', 'plotfile', 'pdf to save to', '/dev/null'))

assocs = io$load(args$assocs)[[args$subset]]

p1 = assocs %>%
    mutate(label = gene_name) %>%
    plt$p_effect("p.value") %>%
    plt$volcano(repel=TRUE, label_top=30) + labs(y="p-value")

highlight = assocs %>%
    head(20) %>%
    mutate(gene_name = b$refactor(gene_name, statistic),
           data = purrr::map(mod, function(m) m$model)) %>%
    select(-mod) %>%
    tidyr::unnest() %>%
    mutate(reads = as.numeric(reads))
smp_aneup = highlight %>%
    arrange(aneup) %>%
    pull(sample) %>%
    unique()
highlight$sample = factor(highlight$sample, levels=smp_aneup)

p21 = ggplot(highlight, aes(x=gene_name, y=sample)) +
    geom_tile(aes(fill=reads), color="white") +
    coord_fixed() +
    viridis::scale_fill_viridis(option="magma", direction=-1) +
    theme(axis.text.x = element_text(size=10, angle=65, hjust=1),
          axis.text.y = element_text(size=10),
          axis.title.x = element_text(size=12),
          legend.position = "left")
p22 = highlight %>%
    select(sample, aneup) %>%
    distinct() %>%
    ggplot(aes(x=aneup, y=sample)) +
        geom_segment(aes(xend=aneup, yend=sample), x=0, color="lightgrey") +
        geom_point() +
        theme(axis.text.x = element_text(size=10),
              axis.title.x = element_text(size=12),
              axis.text.y = element_blank(),
              axis.title.y = element_blank())

pdf(args$plotfile, 10, 8)
print(p1)
print(p21 + p22)
dev.off()
