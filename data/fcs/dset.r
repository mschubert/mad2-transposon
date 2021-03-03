library(Biobase)
library(flowCore)
library(dplyr)
sys = import('sys')

make_gates = function() {
    gates = yaml::read_yaml("fsc-ssc-gates.yaml")
    fsc_ssc = as_tibble(gates$common$debris) * 1e3
#    if (basename(fname) %in% names(gates$sample))
#        fsc_ssc = as_tibble(gates$sample[[basename(fname)]]$debris) * 1e3
    do.call(polygonGate, c(fsc_ssc, list(filterId="debris")))
}

cluster_flowframe = function(ff, transformer, gates, k=1:6, elbow_slope=0.5) {
    ff = flowSet(ff) # can not transform flowFrame?!
    db = flowCore::filter(ff[[1]], gates)
    keep = db@subSet

    # not sure if compensation is already applied, applying multiple times gives no errors
    # ff_comp = flowCore::compensate(ff, spillover(ff)$SPILL)

    options(mc.cores = 1) # flowClust does mclapply, we want to parallelize the whole call
    res = flowWorkspace::GatingSet(ff) %>%
        flowCore::transform(transformer) %>%
        flowWorkspace::gh_pop_get_data() %>% # returns cytoframe
        flowWorkspace::cytoframe_to_flowFrame() %>%
        flowCore::Subset(db) %>%
        flowClust::flowClust(varNames=names(transformer), K=k, B=500)

    if (length(k) > 1) {
        scale01 = function(x) (x - min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE))
        perf = scale01(flowClust::criterion(res, "BIC")) - elbow_slope * scale01(k)
        res = res[[which.max(perf)]]
    }

    # todo: use cluster centers from flowCore and reverse transform?

    df = as_tibble(ff[[1]]@exprs) %>%
        mutate(debris_gate = keep,
               cl = NA)
    df$cl[keep] = res@label
    df$cl = factor(df$cl)
    df
}

sys$run({
    args = sys$cmd$parse(
        opt('d', 'dir', 'fcs file dir', 'FCS files - part 1'),
        opt('o', 'outfile', 'rds', 'FCS files - part 1.rds')
    )

    fcs = list.files(args$dir, pattern="\\.fcs$", recursive=TRUE, full.names=TRUE)
    fs = read.flowSet(fcs)
    bf = flowCore::boundaryFilter(colnames(fs), side="both")
    fs = flowCore::Subset(fs, bf)

    meta = pData(flowCore::parameters(fs[[1]]))
    ffs = flowCore::flowSet_to_list(fs)
    names(ffs) = tools::file_path_sans_ext(names(ffs)) %>%
        sub("Specimen_002_", "", ., fixed=TRUE)

    colors = meta$name[!is.na(meta$desc)]
    biexpTrans = flowWorkspace::flowjo_biexp_trans(channelRange=4096, maxValue=262144, pos=4.5, neg=0, widthBasis=-10)
    trans_colors = flowWorkspace::transformerList(colors, biexpTrans)

    gates = make_gates()

    res = clustermq::Q(cluster_flowframe, ff=ffs,
                       const = list(gates=gates, transformer=trans_colors),
                       n_jobs = 20, memory = 3072, pkgs=c("flowCore"))
    names(res) = names(ffs)

    saveRDS(list(res=res, trans=trans_colors, meta=meta, gates=gates), file=args$outfile)
})
