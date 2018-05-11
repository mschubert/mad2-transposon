library(dplyr)
library(patchwork)
b = import('base')
io = import('io')
plt = import('plot')
sys = import('sys')

plot_volcano = function(assocs) {
    assocs %>%
        mutate(label = gene_name) %>%
        plt$p_effect("p.value") %>%
        plt$volcano(repel=TRUE, label_top=30) + labs(y="p-value")
}

plot_tiles = function(cis, genes=NULL, fill="reads") {
    cis = mutate(cis, sample = b$refactor(sample, aneup))
    aneup = cis %>% select(sample, aneup) %>% distinct()
    if (!is.null(genes))
        cis = mutate(cis, gene_name = factor(gene_name, levels=genes))
    cis = filter(cis, !is.na(aneup) & !is.na(gene_name)) %>%
        tidyr::complete(sample, gene_name, fill=list(reads=0))

    left = ggplot(cis, aes(x=gene_name, y=sample)) +
        geom_tile(aes_string(fill=fill), color="white") +
        coord_fixed() +
        viridis::scale_fill_viridis(option="magma", direction=-1) +
        theme(axis.text.x = element_text(size=10, angle=65, hjust=1),
              axis.text.y = element_text(size=10),
              axis.title.x = element_text(size=12),
              legend.position = "left")
    right = ggplot(aneup, aes(x=aneup, y=sample)) +
        geom_segment(aes(xend=aneup, yend=sample), x=0, color="lightgrey") +
        geom_point() +
        theme(axis.text.x = element_text(size=10),
              axis.title.x = element_text(size=12),
              axis.text.y = element_blank(),
              axis.title.y = element_blank())
    left + right
}

select_highlight = function(assocs) {
    highlight = assocs %>%
        head(20) %>%
        mutate(gene_name = b$refactor(gene_name, statistic),
               data = purrr::map(mod, function(m) m$model)) %>%
        select(-mod) %>%
        tidyr::unnest() %>%
        mutate(reads = as.numeric(reads)) %>%
        group_by(gene_name) %>%
        mutate(rel_reads = reads / max(reads)) %>%
        ungroup()
}

sys$run({
    args = sys$cmd$parse(
        opt('a', 'assocs', 'assocs .RData', 'betareg_adjNreads.RData'),
        opt('s', 'subset', 'assoc object o use', 'hits'),
        opt('p', 'plotfile', 'pdf to save to', '/dev/null'))

    assocs = io$load(args$assocs)[[args$subset]]
    highlight = select_highlight(assocs)

    pdf(args$plotfile, 10, 8)
    print(plot_volcano(assocs))
    print(plot_tiles(highlight, fill="reads"))
    print(plot_tiles(highlight, fill="rel_reads"))
    dev.off()
})
