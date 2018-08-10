library(dplyr)
io = import('io')
sys = import('sys')
gset = import('data/genesets')

test_set = function(set_name) {
    subs = stats$external_gene_name %in% sets[[set_name]]
    len = sum(stats$TTAAs[subs])
    ins = sum(stats$n_smp[subs])
    broom::tidy(poisson.test(ins, len*assocs$n_smp, assocs$ins_rate_genome))
}

args = sys$cmd$parse(
    opt('g', 'gene', 'gene-level poisson', 'poisson.RData'),
    opt('s', 'set', 'gene set file', '{set}.RData'),
    opt('o', 'outfile', 'save results .RData', 'poisson_{set}.RData'))

assocs = io$load(args$gene)
ins = assocs$samples %>%
    group_by(external_gene_name) %>%
    summarize(n_smp = sum(n_ins > 0))
stats = as.data.frame(assocs$gene) %>%
    select(external_gene_name, TTAAs) %>%
    inner_join(ins)

sets = gset$go("mmusculus_gene_ensembl", "external_gene_name", as_list=TRUE) %>%
    gset$filter(min=5, max=200, valid=stats$external_gene_name)

result = sapply(names(sets), test_set, simplify=FALSE) %>%
    dplyr::bind_rows(.id="set") %>%
    select(-(parameter:alternative)) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

save(result, file=args$outfile)
