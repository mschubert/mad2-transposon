library(dplyr)
b = import('base')
io = import('io')
sys = import('sys')
plt = import('plot')
gset = import('data/genesets')

#' Test different insertion statistics of a gene with an external variable
#'
#' This currently does not correct for when a cancer type is both more aneuploid
#' and has more insertions in a given gene.
#'
#' @param dset  data.frame with fields: sample, gene, ext_var[, subset]
#' @param set_name  name of set in sets
#' @param ext_var  differential insertion rate with what, e.g. 'aneuploidy'
#' @param is_type  filter dset with specific 'type' (default: all)
#' @param sets  named list of character vectors for gene sets
#' @return  data.frame with tidy association results
test_set = function(dset, set_name, ext_var, is_type=NA, sets) {
    `%>%` = magrittr::`%>%`
    if (is.na(is_type))
        is_type = unique(dset$type)
    tset = dset %>%
        dplyr::filter(type %in% is_type) %>%
        dplyr::mutate(in_set = external_gene_name %in% sets[[set_name]],
                      ext = !! rlang::sym(ext_var))
    mod = glm(reads ~ external_gene_name + ext + ext:in_set, family=poisson(), data=tset)
    broom::tidy(mod) %>%
        dplyr::mutate(size = sum(tset$reads, na.rm=TRUE)) %>%
        dplyr::filter(term == "ext:in_setTRUE") %>%
        dplyr::select(-term)
}

#' Volcano plot
#'
#' @param res  association result data.frame: external_gene_name, p.value, estimate
#' @return  ggplot2 object
plot_volcano = function(res) {
    res %>%
        mutate(label = set_name) %>%
        plt$p_effect("adj.p", thresh=0.1) %>%
        plt$volcano(p=0.1, label_top=30, base.size=0.2, text.size=2, repel=TRUE)
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'meta+aneuploidy', '../ploidy_compare/analysis_set.RData'),
        opt('c', 'cis_gene', 'poisson cis', '../cis_analysis/poisson.RData'),
        opt('g', 'ext_gene', 'cis assocs RData', '../cis_analysis/ext_gene.RData'), # ignored
        opt('s', 'cis_set', 'cis for sets', 'poisson_set.RData'), # ignored
        opt('o', 'outfile', 'aneuploidy assocs', 'ext_set.RData'),
        opt('p', 'plotfile', 'pdf', 'ext_set.pdf'))

    meta = io$load(args$meta) %>%
        mutate(aneuploidy = pmin(aneuploidy, 0.2))
    dset = io$load(args$cis_gene)
    cis_genes = dset$result %>%
        filter(adj.p < 0.05) %>%
        pull(external_gene_name)
    go = gset$go('mmusculus_gene_ensembl', 'external_gene_name', as_list=TRUE) %>%
        gset$filter(min=3, max=200, valid=cis_genes)
    # valid sets by CIS genes, but then include all genes w/ at least 1 ins?

    # could fit sample rate, gene rate separately and correct with 1 var each?

    aset = dset$samples %>%
        filter(external_gene_name %in% cis_genes) %>%
        select(-n_ins) %>%
        mutate(reads=TRUE) %>%
        tidyr::complete(sample, external_gene_name, fill=list(reads=FALSE)) %>%
        inner_join(meta, by="sample")

    # 1800 sets, 2 hours per set -> 300 jobs = 12 hours (aset @ 160k rows)
    fields = c("aneuploidy", "T-cell", "Myeloid", "Other")
    result = expand.grid(set_name=names(go), ext=fields,
            type=unique(aset$type), stringsAsFactors=FALSE) %>%
        filter(is.na(type) | ext == "aneuploidy") %>%
        mutate(res = clustermq::Q(test_set, const=list(dset=aset, sets=go),
            set_name=set_name, ext_var=ext, is_type=type,
            job_size=5, n_jobs=300, memory=10240)) %>%
        tidyr::unnest() %>%
        mutate(cohens_d = statistic / sqrt(size))

    pdf(args$plotfile)
    for (rn in names(results))
        print(plot_volcano(results[[rn]]) + ggtitle(rn))
    dev.off()

    save(results, file=args$outfile)
})
