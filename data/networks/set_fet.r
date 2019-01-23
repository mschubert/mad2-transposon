library(dplyr)
io = import('io')
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
        opt('n', 'network', 'aracne RData', 'E-GEOD-13159.RData'),
        opt('s', 'setfile', 'gene set RData', '../genesets/human/GO_Biological_Process_2018.RData'),
        opt('o', 'outfile', 'result RData', 'set_fet/GO_Biological_Process_2018.RData')
    )

    net = io$load(args$network)
    sets = io$load(args$setfile) %>%
        gset$filter(min=5, max=500, valid=c(net$Regulator, net$Target))

    result = expand.grid(Regulator = unique(net$Regulator),
                         set = names(sets),
                         stringsAsFactors = FALSE) %>%
        mutate(fet = purrr::map2(Regulator, set, do_fet)) %>%
        tidyr::unnest() %>%
        select(-method, -conf.low, -conf.high, -alternative) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(p.value)

    save(result, file=args$outfile)
})
