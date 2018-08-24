library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
io = import('io')
sys = import('sys')
plt = import('plot')
idmap = import('process/idmap')
gset = import('data/genesets')

plot_gset = function(res, sets, highlight=NULL) {
    test_one = function(set_name) {
        fdata = mutate(cur, in_set = ensembl_gene_id %in% sets[[set_name]])
        mod = try(lm(stat ~ in_set, data=fdata))
        if (class(mod) == "try-error")
            return()
        broom::tidy(mod) %>%
            filter(term == "in_setTRUE") %>%
            select(-term) %>%
            mutate(size = sum(fdata$in_set, na.rm=TRUE))
    }
    cur = res %>% mutate(stat = log2FoldChange / lfcSE)
    result = sapply(names(sets), test_one, simplify=FALSE) %>%
        dplyr::bind_rows(.id="label") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
    result %>%
        plt$p_effect(thresh=0.1) %>%
        plt$volcano(p=0.1, base.size=0.5, label_top=30, repel=TRUE, text.size=2)
}

sys$run({
    args = sys$cmd$parse(
        opt('e', 'de', 'diff expr RData', 'de.RData'),
#        opt('o', 'outfile', 'results RData', 'de_sets.RData'),
        opt('p', 'plotfile', 'pdf', 'de_sets.pdf'))

    res = io$load(args$de)$res
    go = gset$go('mmusculus_gene_ensembl', 'ensembl_gene_id', as_list=TRUE) %>%
        gset$filter(min=5, max=200, valid=rownames(counts))

    pdf(args$plotfile)
    for (name in names(res)) {
        message(name)
        print(plot_gset(res[[name]], go) + ggtitle(name))
    }
    dev.off()

#    save(eset, pca, cis, res, file=args$outfile)
})
