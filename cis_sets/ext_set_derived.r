library(dplyr)
b = import('base')
io = import('io')
sys = import('sys')
plt = import('plot')
gset = import('data/genesets')

#' Test if poisson insertion statistic is different for a set vs rest of genes
#'
#' @param ext_gene  association table from poisson regression
#' @param sets  List of character vectors for gene sets
#' @return  data.frame with associations for all sets
test_sets = function(ext_gene, sets) {
    test_one = function(set_name) {
        cur = mutate(ext_gene, in_set = external_gene_name %in% sets[[set_name]])
        lm(statistic ~ in_set, data=cur) %>%
            broom::tidy() %>%
            filter(term == "in_setTRUE") %>%
            select(-term) %>%
            mutate(size = sum(cur$in_set, na.rm=TRUE))
    }
    result = sapply(names(sets), test_one, simplify=FALSE) %>%
        dplyr::bind_rows(.id="set_name") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
}

#' Volcano plot
#'
#' @param res  association result data.frame: external_gene_name, p.value, estimate
#' @return  ggplot2 object
plot_volcano = function(res) {
    res %>%
        mutate(label = set_name) %>%
        plt$p_effect("adj.p", thresh=0.1) %>%
        plt$volcano(p=0.1, label_top=30, text.size=2, repel=TRUE)
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'meta+aneuploidy', '../ploidy_compare/analysis_set.RData'),
        opt('c', 'cis_gene', 'poisson cis', '../cis_analysis/poisson.RData'), # ignored
        opt('g', 'ext_gene', 'cis assocs RData', '../cis_analysis/ext_gene.RData'),
        opt('s', 'cis_set', 'cis for sets', 'poisson_set.RData'),
        opt('o', 'outfile', 'aneuploidy assocs', 'ext_set_derived.RData'),
        opt('p', 'plotfile', 'pdf', 'ext_set_derived.pdf'))

    meta = io$load(args$meta) %>%
        mutate(aneuploidy = pmin(aneuploidy, 0.2))
    cis_sets = io$load(args$cis_set)$result %>%
        filter(adj.p < 0.05) %>% pull(set)
    dset = io$load(args$ext_gene)
    go = gset$go('mmusculus_gene_ensembl', 'external_gene_name', as_list=TRUE)
    go = go[cis_sets]

    results = lapply(dset, test_sets, sets=go)

    pdf(args$plotfile)
    for (rn in names(results))
        print(plot_volcano(results[[rn]]) + ggtitle(rn))
    dev.off()

    save(results, file=args$outfile)
})