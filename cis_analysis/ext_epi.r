library(dplyr)
b = import('base')
io = import('io')
sys = import('sys')
plt = import('plot')

#' Test different insertion statistics of a gene with an external variable
#'
#' This currently does not correct for when a cancer type is both more aneuploid
#' and has more insertions in a given gene.
#'
#' @param dset  data.frame with fields: sample, epi, ext_var[, subset]
#' @param epi  id to test
#' @param ext_var  differential insertion rate with what, e.g. 'aneuploidy'
#' @param is_type  filter dset with specific 'type' (default: all)
#' @return  data.frame with tidy association results
test_epi = function(dset, epi, ext_var, is_type=NA) {
    `%>%` = magrittr::`%>%`
    if (is.na(is_type))
        is_type = unique(dset$type)
    tset = dset %>%
        dplyr::filter(type %in% is_type, id == epi) %>%
        dplyr::mutate(ext = !! rlang::sym(ext_var))
    broom::tidy(glm(reads ~ ext, family=poisson(), data=tset)) %>%
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
        filter(size >= 2) %>%
        mutate(label = id) %>%
        plt$p_effect("p.value", thresh=0.1) %>%
        plt$volcano(p=0.1, label_top=30, repel=TRUE) +
            ylab("non-adj. p-value")
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'meta+aneuploidy', '../ploidy_compare/analysis_set.RData'),
        opt('p', 'poisson', 'cis assocs RData', 'poisson_epi.RData'),
        opt('o', 'outfile', 'external assocs RData', 'ext_epi.RData'),
        opt('p', 'plotfile', 'pdf', 'ext_epi.pdf'))

    meta = io$load(args$meta) %>%
        mutate(aneuploidy = pmin(aneuploidy, 0.2))
    dset = io$load(args$poisson)
    ids = with(dset, intersect(samples$id, result$id))

    aset = dset$samples %>%
        select(-n_ins) %>%
        filter(id %in% ids) %>%
        mutate(reads=TRUE) %>%
        tidyr::complete(sample, id, fill=list(reads=FALSE)) %>%
        inner_join(dset$sample_rates %>% select(sample, n_ins, n_reads, rate)) %>%
        inner_join(meta %>% select(sample, aneuploidy, type,
                                   `T-cell`, Myeloid, Other))

    fields = c("aneuploidy", "T-cell", "Myeloid", "Other")
    result = expand.grid(id=ids, ext=fields,
            type=unique(aset$type), stringsAsFactors=FALSE) %>%
        filter(is.na(type) | ext == "aneuploidy") %>%
        mutate(res = clustermq::Q(test_epi, const=list(dset=aset),
            epi=id, ext_var=ext, is_type=type,
            job_size=50, n_jobs=5, memory=1024)) %>%
        tidyr::unnest() %>%
        mutate(cohens_d = statistic / sqrt(size))

    results = result %>%
        mutate(subset = ifelse(is.na(type), ext, paste0(type, ":", ext))) %>%
        select(-ext, -type) %>%
        split(.$subset)

    pdf(args$plotfile)
    for (rn in names(results)) {
        message(rn)
        print(plot_volcano(results[[rn]]) + ggtitle(rn))
    }
    dev.off()

    save(results, file=args$outfile)
})
