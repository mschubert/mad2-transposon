library(flowCore)
library(dplyr)
sys = import('sys')

process_one = function(fname, cluster=TRUE) {
    message(fname)
#    ff = read.FCS(fname)
    ff = read.flowSet(fname)
    bn = tools::file_path_sans_ext(basename(fname))
    bn = sub("Specimen_002_", "", bn, fixed=TRUE)

    gates = yaml::read_yaml("fsc-ssc-gates.yaml")
    fsc_ssc = as_tibble(gates$common$debris) * 1e3
    if (basename(fname) %in% names(gates$sample))
        fsc_ssc = as_tibble(gates$sample[[basename(fname)]]$debris) * 1e3
    debris = do.call(polygonGate, c(fsc_ssc, list(filterId="debris")))

    meta = flowCore::parameters(ff[[1]]) %>% Biobase::pData()
    db = flowCore::filter(ff[[1]], debris)
    keep = db@subSet

#    df = as_tibble(Subset(ff, db)@exprs) %>%
#        dplyr::filter(`SSC-A` < 2.5e5) # todo: is this caught by singlets [FSC-H FSC-A] gate?
#    df = df[rowSums(df < 0) == 0,]

    fields = c(na.omit(meta$desc))
    colors = meta$name[!is.na(meta$desc)]

    # not sure if compensation is already applied, applying multiple times gives no errors
    # ff_comp = flowCore::compensate(ff, spillover(ff)$SPILL)

    biexpTrans = flowWorkspace::flowjo_biexp_trans(channelRange=4096, maxValue=262144, pos=4.5, neg=0, widthBasis=-10)
    trans_colors = flowWorkspace::transformerList(colors, biexpTrans)

    options(mc.cores = 6)
    res = flowWorkspace::GatingSet(ff) %>%
        flowCore::transform(trans_colors) %>%
        flowWorkspace::gh_pop_get_data() %>% # returns cytoframe
        flowWorkspace::cytoframe_to_flowFrame() %>%
        flowCore::Subset(db) %>%
        flowClust::flowClust(varNames=colors, K=1:6)
#    res[[1]]@mu # this is transformed space
    best = which.max(scale(flowClust::criterion(res, "BIC")) - scale(1)) # elbow method

    df = as_tibble(ff[[1]]@exprs)
    colnames(df) = ifelse(is.na(meta$desc), meta$name, meta$desc)
    df = df %>%
        mutate(debris_gate = keep,
               cl = NA)
    df$cl[keep] = factor(res[[best]]@label)
    df
}

sys$run({
    args = sys$cmd$parse(
        opt('d', 'dir', 'fcs file dir', 'FCS files - part 1'),
        opt('o', 'outfile', 'rds', 'FCS files - part 1.rds')
    )

    set.seed(120587)
    fcs = list.files(args$dir, pattern="\\.fcs$", recursive=TRUE, full.names=TRUE)

    # fcs = sample(fcs, 3)

    res = lapply(fcs, function(x) try(plot_one(x))) %>%
        setNames(fcs)
    errs = sapply(res, class) == "try-error"
    if (any(errs)) {
        warning("Samples failed: ", paste(names(res)[errs], collapse=", "))
        res = res[!errs]
    }

    saveRDS(res, file=args$outfile)
})
