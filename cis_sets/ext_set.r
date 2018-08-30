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
    n_ins = sum(with(tset, reads[in_set]))
    if (n_ins < 5)
        return()
    mod6 = lme4::glmer(reads ~ ext + ext:in_set + (1|sample) + 
        (ext|type) + (1|external_gene_name), family=poisson(), data=tset,
        control=lme4::glmerControl(optCtrl=list(method="L-BFGS-B")))
    broom::tidy(mod) %>%
        dplyr::filter(term == "ext:in_setTRUE") %>%
        dplyr::mutate(n_genes = length(sets[[set_name]]),
                      n_ins = n_ins) %>%
        dplyr::select(n_genes, n_ins, estimate:p.value)
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
        opt('f', 'sets', 'RData for gene set', '../data/genesets/KEA_2015.RData'),
        opt('o', 'outfile', 'aneuploidy assocs', 'ext_set/KEA_2015.RData'),
        opt('p', 'plotfile', 'pdf', 'ext_set/KEA_2015.pdf'))

    meta = io$load(args$meta) %>%
        mutate(aneuploidy = pmin(aneuploidy, 0.2))
    dset = io$load(args$cis_gene)
    sets = io$load(args$sets)
        gset$filter(min=5, valid=dset$result$external_gene_name)

    aset = dset$samples %>%
        select(-n_ins) %>%
        mutate(reads=TRUE) %>%
        tidyr::complete(sample, external_gene_name, fill=list(reads=FALSE)) %>%
        inner_join(meta, by="sample")

    # 100 @ 20 jobs = 30 min; 25k @ 200 (x25) -> 12 hours
    fields = c("aneuploidy", "T-cell", "Myeloid", "Other")
    result = expand.grid(set_name=names(sets), ext=fields,
            type=unique(aset$type), stringsAsFactors=FALSE) %>%
        filter(is.na(type) | ext == "aneuploidy") %>%
        mutate(res = clustermq::Q(test_set, const=list(dset=aset, sets=sets),
            set_name=set_name, ext_var=ext, is_type=type,
            job_size=5, n_jobs=300, memory=3072)) %>%
        tidyr::unnest()

    results = res %>%
        mutate(label = ifelse(is.na(type), ext, paste(type, ext, sep=":")),
               size = n_ins,
               adj.p = p.adjust(p.value, method="fdr"))
    results = split(results, results$label)

    pdf(args$plotfile)
    for (rn in names(results))
        print(plot_volcano(results[[rn]]) + ggtitle(rn))
    dev.off()

    save(result, file=args$outfile)
})
