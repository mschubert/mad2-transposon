library(dplyr)
b = import('base')
io = import('io')
sys = import('sys')
enr = import('tools/enrichr')
msd = import('tools/msigdb')

#' Score gene expression using KS statistic on different sets
#'
#' @param data  A data.frame with fields: 'Gene'
#' @param sets  A list of character vectors
#' @return     
score = function(data, sets) {
    # KS-statistic p-value per set (and mean of change) [or piano]
    # using enrichr$run is quite unstable with different cutoffs
    do_ks = function(set) {
        matches = which(data$Gene %in% set)
        data.frame(effect = data %>% filter(Gene %in% set) %$% mean(residual),
                   p.value = ks.test(matches, (1:nrow(data))[-matches])$p.value,
                   size = length(set))
    }
    scores = lapply(sets, function(s) do_ks(s) %catch% NULL) %>%
        b$omit$null() %>%
        df$add_name_col("set", bind=TRUE) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
}

sys$run({
    args = sys$cmd$parse(
        opt('g', 'geneset', 'Identifier of the gene set'),
        opt('o', 'outfile', 'File to save gene set to'))

    if (args$geneset %in% enr$dbs()$name)
        sets = enr$genes(args$geneset)
    else if (args$geneset %in% msd$dbs())
        sets = msd$genes(args$geneset)
    else
        stop("invalid gene set: ", args$geneset)

    io$save(sets, file=args$outfile)
})
