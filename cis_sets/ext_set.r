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
test_set = function(dset, set_name, ext_var, is_type=NA, sets, min_n=5) {
    library(dplyr)
    if (!is.na(is_type))
        dset = dset %>% filter(type %in% is_type)

    # regressing out ins_total yields (almost) no signif assocs
    # i.e. diff gene sets are mainly driven by overall number
    if (ext_var == "aneuploidy" && is.na(is_type))
        fml = ins_set ~ type + ins_total + ext
    else
        fml = ins_set ~ ins_total + ext

    tset = dset %>%
        mutate(in_set = external_gene_name %in% sets[[set_name]],
               ext = !! rlang::sym(ext_var)) %>%
        group_by(sample) %>%
        summarize(type = unique(type),
                  ins_set = sum(reads[in_set]),
                  ins_total = sum(reads),
                  ext = unique(ext))

    n_genes = length(sets[[set_name]])
    n_ins = sum(tset$ins_set)
    if (n_genes < min_n || n_ins < min_n)
        return(data.frame(n_genes = n_genes, n_ins=n_ins))

    mod = glm(fml, family=poisson(), data=tset)
    broom::tidy(mod) %>%
        filter(term == "ext") %>%
        mutate(n_genes = n_genes, n_ins = n_ins) %>%
        select(n_genes, n_ins, estimate:p.value)
}

#' Volcano plot
#'
#' @param res  association result data.frame: external_gene_name, p.value, estimate
#' @return  ggplot2 object
plot_volcano = function(res) {
    res %>%
        mutate(label = set_name, size = n_genes) %>%
#        plt$p_effect("adj.p", thresh=0.1) %>%
        plt$p_effect("p.value", thresh=0.05) %>%
        plt$volcano(p=0.05, label_top=30, base.size=0.1, text.size=2, repel=TRUE) +
            ylab("p-value")
}

sys$run({
    args = sys$cmd$parse(
        opt('m', 'meta', 'meta+aneuploidy', '../ploidy_compare/analysis_set.RData'),
        opt('c', 'cis_gene', 'poisson cis', '../cis_analysis/poisson.RData'),
        opt('g', 'ext_gene', 'cis assocs RData', '../cis_analysis/ext_gene.RData'), # ignored
        opt('s', 'cis_set', 'cis for sets', 'poisson_set/GO_Biological_Process_2018.RData'), # ignored
        opt('f', 'sets', 'genes per set', '../data/genesets/mouse/GO_Biological_Process_2018.RData'),
        opt('o', 'outfile', 'aneuploidy assocs', 'ext_set/GO_Biological_Process_2018.RData'),
        opt('p', 'plotfile', 'pdf', 'ext_set/GO_Biological_Process_2018.pdf'))

    meta = io$load(args$meta) %>%
        mutate(aneuploidy = pmin(aneuploidy, 0.2))
    dset = io$load(args$cis_gene)
    sets = io$load(args$sets) %>%
        gset$filter(min=5, valid=dset$result$external_gene_name)
#    cis_sets = io$load(args$cis_set)$result %>%
#        filter(abs(estimate) > 2 & adj.p < 0.01) %>%
#        pull(set)
#    sets = sets[cis_sets]

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
            job_size=500, n_jobs=100, memory=1024, chunk_size=1)) %>%
        tidyr::unnest() %>%
        filter(!is.na(p.value)) %>%
        mutate(label = ifelse(is.na(type), ext, paste(type, ext, sep=":"))) %>%
        group_by(label) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        ungroup() %>%
        arrange(adj.p)

    results = split(result, result$label)
    pdf(args$plotfile)
    for (rn in names(results))
        print(plot_volcano(results[[rn]]) + ggtitle(rn))
    dev.off()

    save(result, file=args$outfile)
})
