library(dplyr)
library(patchwork)
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
    plt$p_effect("adj.p") %>%
    plt$volcano(repel=TRUE, label_top=30)

highlight = assocs %>%
    head(20) %>%
    mutate(gene_name = factor(gene_name, levels=gene_name),
           data = purrr::map(mod, function(m) m$model)) %>%
    select(-mod) %>%
    tidyr::unnest()
smp_aneup = highlight %>%
    arrange(aneup) %>%
    pull(sample) %>%
    unique()
highlight$sample = factor(highlight$sample, levels=smp_aneup)

p21 = ggplot(highlight, aes(x=gene_name, y=sample)) +
    geom_tile(aes(fill=reads)) +
    coord_fixed() +
    theme(axis.text.x = element_text(angle=65, hjust=1),
          legend.position = "left")
p22 = highlight %>%
    select(sample, aneup) %>%
    distinct() %>%
    ggplot(aes(x=aneup, y=sample)) +
        geom_segment(aes(xend=aneup, yend=sample), x=0, color="lightgrey") +
        geom_point() +
        theme(axis.text.y = element_blank(),
              axis.title.y = element_blank())

pdf(args$plotfile, 10, 8)
print(p1)
print(p21 + p22)
dev.off()
