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
    broom::tidy(poisson.test(ins, len*assocs$n_smp, assocs$ins_rate_genome)) %>%
        mutate(size = sum(subs),
               n_smps = length(unique(smps)),
               samples = list(sample=smps))
}

args = sys$cmd$parse(
    opt('i', 'infile', 'gene-level poisson', 'poisson.RData'),
#    opt('s', 'set', 'gene set file', '{set}.RData'),
    opt('o', 'outfile', 'save results .RData', 'poisson_set.RData'),
    opt('p', 'plotfile', 'pdf', 'poisson_set.pdf'))

assocs = io$load(args$infile)
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

sets = gset$go("mmusculus_gene_ensembl", "external_gene_name", as_list=TRUE) %>%
    gset$filter(min=5, max=200, valid=stats$external_gene_name)

result = sapply(names(sets), test_set, simplify=FALSE) %>%
    dplyr::bind_rows(.id="set") %>%
    select(-(parameter:alternative)) %>%
    select(set, size, n_smps, everything()) %>%
    mutate(estimate = (estimate - assocs$ins_rate_genome) / assocs$ins_rate_genome,
           adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

samples = result %>%
    select(set, samples) %>%
    tidyr::unnest() %>%
    group_by(samples, set) %>%
    summarize(n_ins = n()) %>%
    ungroup() %>%
    dplyr::rename(sample = samples)

result = result %>%
    select(-samples)

p = result %>%
    mutate(label = ifelse((40/6)*estimate-40 > log10(adj.p) | estimate < -0.5 ,
                          set, NA)) %>%
    plt$p_effect() %>%
    plt$volcano(text.size=2, label_top=Inf, repel=TRUE)

pdf(args$plotfile, 15, 6)
print(p)
dev.off()

save(result, samples, file=args$outfile)
