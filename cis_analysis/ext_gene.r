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
#' @param dset  data.frame with fields: sample, gene, ext_var[, subset]
#' @param gene  external_gene_name to test
#' @param ext_var  differential insertion rate with what, e.g. 'aneuploidy'
#' @param is_type  filter dset with specific 'type' (default: all)
#' @return  data.frame with tidy association results
test_gene = function(dset, gene, ext_var, is_type=NA) {
    `%>%` = magrittr::`%>%`
    if (is.na(is_type)) {
        is_type = unique(dset$type)
        fml = reads ~ type + ext
    } else
        fml = reads ~ ext
    tset = dset %>%
        dplyr::filter(type %in% is_type, external_gene_name == gene) %>%
        dplyr::mutate(ext = !! rlang::sym(ext_var))
    broom::tidy(glm(fml, family=poisson(), data=tset)) %>%
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
        mutate(label = external_gene_name) %>%
        plt$p_effect("p.value", thresh=0.1) %>%
        plt$volcano(p=0.1, label_top=30, repel=TRUE) +
            ylab("non-adj. p-value")
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'meta+aneuploidy', '../ploidy_compare/analysis_set.RData'),
        opt('p', 'poisson', 'cis assocs RData', 'poisson.RData'),
        opt('o', 'outfile', 'external assocs RData', 'ext_gene.RData'),
        opt('p', 'plotfile', 'pdf', 'ext_gene.pdf'))

    meta = io$load(args$meta)
    dset = io$load(args$poisson)
    genes = with(dset, intersect(samples$external_gene_name,
                                 result$external_gene_name))

    aset = dset$samples %>%
        select(-n_ins) %>%
        filter(external_gene_name %in% genes) %>%
        mutate(reads=TRUE) %>%
        tidyr::complete(sample, external_gene_name, fill=list(reads=FALSE)) %>%
        inner_join(dset$sample_rates %>% select(sample, n_ins, n_reads, rate)) %>%
        inner_join(meta %>% select(sample, aneuploidy, type,
                                   `T-cell`, Myeloid, Other))

    fields = c("aneuploidy", "T-cell", "Myeloid", "Other")
    result = expand.grid(external_gene_name=genes, ext=fields,
            type=unique(aset$type), stringsAsFactors=FALSE) %>%
        filter(is.na(type) | ext == "aneuploidy") %>%
        mutate(res = clustermq::Q(test_gene, const=list(dset=aset),
            gene=external_gene_name, ext_var=ext, is_type=type,
            job_size=50, n_jobs=20, memory=1024, fail_on_error=FALSE)) %>%
        filter(sapply(res, class) != "error") %>%
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
