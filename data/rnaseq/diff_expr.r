library(lme4)
b = import('base')
io = import('io')
st = import('stats')
gs = import('data/genesets')

dset = io$load('assemble.RData')
expr = dset$expr 
expr = expr[narray::map(expr, along=2, sd) > 0.15,]
expr = narray::map(expr, along=2, scale)

go = gs$go(dset="mmusculus_gene_ensembl", genes="ensembl_gene_id") %>%
    filter(ensembl_gene_id %in% rownames(expr)) %>%
    transmute(gene = ensembl_gene_id,
              set = paste(id, name)) %>%
    unstack()
size = sapply(go, length)
go = go[size >= 5 & size <= 200]

#TODO: make constrasts matrix (0, 1, NA)

condition = dset$idx$tissue=="spleen"

test_sets = function(sets, condition) {
    test_set = function(set) {
        t.test(expr[set, condition], expr[set, !condition]) %>%
            broom::tidy() %>%
            mutate(size = length(set)) %>%
            select(size, estimate, statistic, p.value)
    }
    stopifnot(is.logical(condition))
    lapply(sets, function(s) test_set(s) %catch% NA) %>%
        bind_rows() %>%
        mutate(set = names(sets),
               adj.p = p.adjust(p.value, method="fdr")) %>%
        select(set, everything())
}

#TODO: test condition (each contrast), all sets
