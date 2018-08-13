library(dplyr)
b = import('base')
io = import('io')
sys = import('sys')
plt = import('plot')

#' source sample IDs are both 123S and T567
fix_ids = function(x) {
    letter = tolower(sub("[0-9]+", "", x))
    nums = sub("[ST]", "", x)
    paste0(nums, letter)
}

#' KS statistic
ks_test = function(data, field) {
    mod = ks.test(data[[field]][data$reads], data[[field]][!data$reads])
    mod$model = select(data, sample, !!rlang::sym(field), reads)
    mod %>%
        broom::tidy() %>%
        mutate(estimate = mean(data[[field]][data$reads]) - mean(data[[field]]),
               statistic = sign(estimate) / p.value,
               size = sum(data$reads),
               mod = list(mod))
}

#' Poisson regression
poisson_reg = function(data, field) {
    fml = formula(paste("reads ~", field))
    mod = glm(fml, family=poisson(), data=data)
    mod$model = select(data, sample, !!rlang::sym(field), reads)
    mod %>%
        broom::tidy() %>%
        filter(term == field) %>%
        select(-term) %>%
        mutate(size = sum(data$reads),
               mod = list(mod))
}

#' Run given test for a field
do_test = function(dset, test, field) {
    test_fun = switch(args$test, 'ks'=ks_test, 'poisson'=poisson_reg,
                      stop("invalid 'test' argument"))

    dset %>%
        group_by(external_gene_name) %>%
        tidyr::nest() %>%
        mutate(fit = purrr::map(data, test_fun, field=field)) %>%
        select(-data) %>%
        tidyr::unnest() %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        select(external_gene_name, size, estimate, statistic, p.value, adj.p) %>%
        arrange(adj.p, p.value)
}

sys$run({
    args = sys$cmd$parse(
        opt('a', 'aneup', 'aneuploidy .tsv', '../ploidy_compare/analysis_set.RData'),
        opt('p', 'poisson', 'cis assocs RData', 'poisson.RData'),
        opt('t', 'test', 'ks/poisson', 'ks'),
        opt('o', 'outfile', 'aneuploidy assocs', 'aneup_assocs/ks.RData'),
        opt('p', 'plotfile', 'pdf', 'aneup_assocs/ks.pdf'))

    aneup = io$load(args$aneup)
    dset = io$load(args$poisson)
    genes = dset$result %>%
        filter(adj.p < 1e-5) %>%
        pull(external_gene_name)

    aset = dset$samples %>%
        select(-n_ins) %>%
        filter(external_gene_name %in% genes) %>%
        mutate(reads=TRUE) %>%
        tidyr::complete(sample, external_gene_name, fill=list(reads=FALSE)) %>%
        inner_join(aneup, by="sample") %>%
        mutate(type_Tcell = ifelse(type == "T-cell", 1, 0),
               type_Myeloid = ifelse(type == "Myeloid", 1, 0),
               type_Other = ifelse(type == "Other", 1, 0),
               aneup_Tcell = ifelse(type == "T-cell", aneup, 0),
               aneup_Myeloid = ifelse(type == "Myeloid", aneup, 0),
               aneup_Other = ifelse(type == "Other", aneup, 0))

    plot_volcano = . %>%
        mutate(label = external_gene_name) %>%
        plt$p_effect(thresh=0.1) %>%
        plt$volcano(p=0.1, label_top=30, repel=TRUE)

#    result = sapply()
#    plots = lapply(result, plot_volcano)

    pdf(args$plotfile)
    print(plot_volcano(do_test(aset, args$test, "aneup")) + ggtitle("aneup"))
    print(plot_volcano(do_test(aset, args$test, "type_Tcell")) + ggtitle("T-cell"))
    print(plot_volcano(do_test(aset, args$test, "type_Myeloid")) + ggtitle("Myeloid"))
    print(plot_volcano(do_test(aset, args$test, "type_Other")) + ggtitle("Other"))
    print(plot_volcano(do_test(aset, args$test, "aneup_Tcell")) + ggtitle("aneup T-cell"))
    print(plot_volcano(do_test(aset, args$test, "aneup_Myeloid")) + ggtitle("aneup Myeloid"))
    print(plot_volcano(do_test(aset, args$test, "aneup_Other")) + ggtitle("aneup Other"))
    dev.off()

    save(genes, file=args$outfile) # rather save result
})
