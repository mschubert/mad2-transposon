sys = import('sys')
gdsc = import('data/gdsc')

load_set = function(yset, hl) {
    x = readRDS(file.path("../data/geneset_scores/gdsc/", paste0(yset, ".rds")))
    re = x[hl[[yset]], ,drop=FALSE]
}

create_edf = function(dset, hl, sets, genes) {
    tissues = gdsc$tissues(dset$tissues)
    expr = gdsc$basal_expression()
    expr = expr[intersect(rownames(expr), genes),
                intersect(colnames(expr), names(tissues))]
    eset = sets[,intersect(colnames(expr), names(tissues))]
    #eset = narray::map(eset, along=2, scale)
    expr = narray::stack(list(expr, eset), along=1)
    colnames(expr) = gdsc$cosmic$id2name(colnames(expr))
    aneup = gdsc$aneuploidy() %>% select(cline=cell_line, aneuploidy)

    edf = reshape2::melt(expr) %>%
        transmute(gene=Var1, cline=as.character(Var2), expr=value,
                  is_gene=ifelse(gene %in% genes, "gene", "geneset")) %>%
        mutate(label = ifelse(cline %in% dset$clines, cline, NA)) %>%
        left_join(aneup)
}

plot_edf = function(edf) {
    p = ggplot(edf, aes(x=gene, y=expr)) +
        ggbeeswarm::geom_quasirandom(aes(alpha=is.na(label), size=aneuploidy)) +
        ggrepel::geom_text_repel(aes(label=label)) +
        guides(alpha=FALSE) +
        facet_wrap(~ is_gene, ncol=1, scales="free") +
        scale_alpha_manual(values=c(1, 0.1)) +
        theme(axis.text.x = element_text(angle=10, hjust=1))
}

sys$run({
    args = sys$cmd$parse(
        opt('h', 'highlight', 'yaml', 'highlight.yaml'),
        opt('d', 'dset', 'gdsc|ccle.yaml', 'gdsc.yaml'),
        opt('p', 'plotfile', 'pdf', 'gdsc.pdf'))

    hl = yaml::read_yaml(args$highlight)
    dset = yaml::read_yaml(args$dset)
    sets = setdiff(names(hl), "genes") %>%
        lapply(load_set, hl=hl) %>%
        narray::stack(along=2)

    plots = lapply(dset, function(x) plot_edf(create_edf(x, hl, sets, hl$genes)))

    pdf(args$plotfile, 14, 12)
    for (i in seq_along(plots))
        print(plots[[i]] + ggtitle(names(dset)[i]))
    dev.off()
})
