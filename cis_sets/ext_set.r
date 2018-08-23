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
        dplyr::filter(type %in% is_type,
                      external_gene_name %in% sets[[set_name]]) %>%
        dplyr::mutate(ext = !! rlang::sym(ext_var))
    mod = glm(reads ~ external_gene_name + ext, family=poisson(), data=tset)
    broom::tidy(mod) %>%
        dplyr::mutate(size = sum(tset$reads, na.rm=TRUE)) %>%
        dplyr::filter(term == "ext") %>%
        dplyr::select(-term)
}

#' Volcano plot
#'
#' @param res  association result data.frame: external_gene_name, p.value, estimate
#' @return  ggplot2 object
plot_volcano = function(res) {
    res %>%
        mutate(label = set_name) %>%
        plt$p_effect(thresh=0.1) %>%
        plt$volcano(p=0.1, label_top=30, base.size=0.2, text.size=2, repel=TRUE) +
            ylab("p-value (non-adj.)")
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'meta+aneuploidy', '../ploidy_compare/analysis_set.RData'),
        opt('p', 'poisson', 'cis assocs RData', '../cis_analysis/poisson.RData'),
        opt('o', 'outfile', 'aneuploidy assocs', 'aneup_assocs.RData'),
        opt('p', 'plotfile', 'pdf', 'aneup_assocs.pdf'))

    meta = io$load(args$meta) %>%
        mutate(aneuploidy = pmin(aneuploidy, 0.2))
    dset = io$load(args$poisson)
    go = gset$go('mmusculus_gene_ensembl', 'external_gene_name', as_list=TRUE) %>%
        gset$filter(min=5, max=200, valid=dset$result$external_gene_name)

    aset = dset$samples %>%
        select(-n_ins) %>%
        mutate(reads=TRUE) %>%
        tidyr::complete(sample, external_gene_name, fill=list(reads=FALSE)) %>%
        inner_join(meta, by="sample")

    fields = c("aneuploidy", "T-cell", "Myeloid", "Other")
    result = expand.grid(set_name=names(go), ext=fields,
            type=unique(aset$type), stringsAsFactors=FALSE) %>%
        filter(is.na(type) | ext == "aneuploidy") %>%
        mutate(res = clustermq::Q(test_set, const=list(dset=aset, sets=go),
            set_name=set_name, ext_var=ext, is_type=type,
            job_size=50, n_jobs=20, memory=5120)) %>%
        tidyr::unnest() %>%
        mutate(cohens_d = statistic / sqrt(size))

    results = result %>%
        mutate(subset = ifelse(is.na(type), ext, paste0(type, ":", ext))) %>%
        select(-ext, -type) %>%
        split(.$subset)

    pdf(args$plotfile)
    for (rn in names(results))
        print(plot_volcano(results[[rn]]) + ggtitle(rn))
    dev.off()

    save(results, file=args$outfile)
})
