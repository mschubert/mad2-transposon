library(dplyr)
io = import('io')
sys = import('sys')
gset = import('data/genesets')
plt = import('plot')

test_set = function(set_name) {
    subs = stats$external_gene_name %in% sets[[set_name]]
    len = sum(stats$TTAAs[subs])
    ins = sum(stats$n_smp[subs])
    smps = unlist(stats$sample[subs])
    if (length(smps) == 0)
        return()
    broom::tidy(poisson.test(ins, len*n_smp, ins_rate_genome)) %>%
        mutate(p.value = ifelse(is.logical(p.value), 1, p.value),
               size = sum(subs),
               n_smps = length(unique(smps)),
               samples = list(sample=smps))
}

args = sys$cmd$parse(
    opt('i', 'infile', 'gene-level poisson', '../cis_analysis/poisson.RData'),
    opt('s', 'sets', 'gene set file', '../data/genesets/KEA_2015.RData'),
    opt('o', 'outfile', 'save results .RData', 'poisson_set.RData'),
    opt('p', 'plotfile', 'pdf', 'poisson_set.pdf'))

assocs = io$load(args$infile)
n_smp = nrow(assocs$sample_rates)
ins_rate_genome = mean(assocs$sample_rates$rate)
cis_genes = assocs$result %>%
    filter(adj.p < 1e-5) %>%
    pull(external_gene_name)
ins = assocs$samples %>%
    group_by(external_gene_name) %>%
    summarize(n_smp = sum(n_ins > 0),
              sample = list(unique(sample))) %>%
    filter(!external_gene_name %in% cis_genes)
stats = as.data.frame(assocs$gene) %>%
    select(external_gene_name, TTAAs) %>%
    inner_join(ins)

sets = io$load(args$sets) %>%
    gset$filter(min=5, max=200)

result = sapply(names(sets), test_set, simplify=FALSE) %>%
    dplyr::bind_rows(.id="set") %>%
    select(-(parameter:alternative)) %>%
    select(set, size, n_smps, everything()) %>%
    mutate(estimate = log2(estimate / ins_rate_genome),
           adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

samples = result %>%
    select(set, samples) %>%
    filter(sapply(samples, is.character) == TRUE) %>%
    tidyr::unnest() %>%
    group_by(samples, set) %>%
    summarize(n_ins = n()) %>%
    ungroup() %>%
    dplyr::rename(sample = samples)

result = result %>%
    select(-samples)

p = result %>% # y = y0 - y0/x0 * x, all *(-1)
#    mutate(label = ifelse((60/3.3)*estimate-60 > log10(adj.p) | estimate < -1, set, NA)) %>%
    mutate(label = set) %>%
    plt$p_effect("adj.p", thresh=0.1) %>%
    plt$volcano(p=0.1, base.size=0.2, text.size=2, label_top=30, repel=TRUE) +
        xlab("log2 fold change poisson rate")

pdf(args$plotfile, 8, 6)
print(p)
dev.off()

save(result, samples, file=args$outfile)
