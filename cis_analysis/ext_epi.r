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

#' Plot distribution of insertion-feature distances vs. external variable
#'
#' @param meta  metadata df
#' @param rn    column name that is external var
#' @param dist  distance measure
#' @return      ggplot2 object
plot_distance_distr = function(meta, rn, dist) {
    df = as_tibble(as.data.frame(dist)) %>%
        inner_join(meta, by="sample") %>%
        dplyr::select(sample, id, distance, !! rlang::sym(rn)) %>%
        mutate(sample = forcats::fct_reorder(sample, !! rlang::sym(rn)),
               distance = log10(distance+1))

    mod = broom::tidy(lm(distance ~ as.integer(sample), data=df)) %>%
        filter(term == "as.integer(sample)")

    ggplot(df, aes(x=sample, y=log10(distance+1))) +
        geom_boxplot(outlier.shape=NA) +
        labs(title = rn,
             subtitle = sprintf("est=%.2f (p=%.2g)", mod$estimate, mod$p.value)) +
        theme(axis.text.x = element_text(angle=90))
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'meta+aneuploidy', '../ploidy_compare/analysis_set.rds'),
        opt('p', 'poisson', 'cis assocs rds', 'poisson_epi/enhancer.rds'),
        opt('o', 'outfile', 'external assocs rds', 'ext_epi/enhancer.rds'),
        opt('p', 'plotfile', 'pdf', 'ext_epi/enhancer.pdf')
    )

    meta = readRDS(args$meta)$meta %>%
        mutate(aneuploidy = pmin(aneuploidy, 0.2))
    dset = readRDS(args$poisson)
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
        if (rn == "aneuploidy")
            print(plot_distance_distr(meta, rn, dset$dist))
    }
    dev.off()

    saveRDS(results, file=args$outfile)
})
