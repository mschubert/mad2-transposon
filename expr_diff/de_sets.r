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
        fdata = mutate(cur, in_set = gene_name %in% sets[[set_name]])
        mod = try(lm(stat ~ in_set, data=fdata))
        if (class(mod) == "try-error")
            return()
        broom::tidy(mod) %>%
            filter(term == "in_setTRUE") %>%
            select(-term) %>%
            mutate(size = sum(fdata$in_set, na.rm=TRUE))
    }
    cur = res %>%
        mutate(stat = log2FoldChange / lfcSE,
               gene_name = idmap$gene(ensembl_gene_id,
                    to="external_gene_name", dset="mmusculus_gene_ensembl"))
    result = sapply(names(sets), test_one, simplify=FALSE) %>%
        dplyr::bind_rows(.id="label") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
    result %>%
        plt$p_effect("adj.p", thresh=0.1) %>%
        plt$volcano(p=0.1, base.size=0.1, label_top=30, repel=TRUE, text.size=2)
}

sys$run({
    args = sys$cmd$parse(
        opt('e', 'de', 'diff expr RData', 'de.RData'),
        opt('s', 'sets', 'gene set RData', '../data/genesets/GO_Biological_Process_2017b.RData'),
        opt('p', 'plotfile', 'pdf', 'de_sets.pdf'))

    res = io$load(args$de)$res
    sets = io$load(args$sets) %>%
        gset$filter(min=5, valid=rownames(counts))

    pdf(args$plotfile)
    for (name in names(res)) {
        message(name)
        print(plot_gset(res[[name]], sets) + ggtitle(name))
    }
    dev.off()
})
