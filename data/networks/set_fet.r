library(dplyr)
sys = import('sys')
gset = import('data/genesets')

do_fet = function(rname, sname) {
    rset = net %>% filter(Regulator == rname) %>% pull(Target) %>% unique()
    rtotal = length(setdiff(net$Target, rset))
    nset = sets[[sname]]
    ng = length(setdiff(unlist(sets), nset))
    fisher.test(matrix(c(length(nset), ng, length(rset), rtotal), ncol=2)) %>%
        broom::tidy()
}

sys$run({
    args = sys$cmd$parse(
        opt('n', 'network', 'rds', 'aracne_E-GEOD-13159.rds'),
        opt('s', 'setfile', 'rds', '../genesets/human/GO_Biological_Process_2020.rds'),
        opt('o', 'outfile', 'rds', 'set_fet/GO_Biological_Process_2020.rds')
    )

    net = readRDS(args$network)
    sets = readRDS(args$setfile) %>%
        gset$filter(min=5, max=500, valid=c(net$Regulator, net$Target))

    result = expand.grid(Regulator = unique(net$Regulator),
                         set = names(sets),
                         stringsAsFactors = FALSE) %>%
        mutate(fet = purrr::map2(Regulator, set, do_fet)) %>%
        tidyr::unnest() %>%
        select(-method, -conf.low, -conf.high, -alternative) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(p.value)

    saveRDS(result, file=args$outfile)
})
