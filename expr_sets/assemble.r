library(dplyr)
library(cowplot)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('d', 'dset', 'expr RData', '../expr_diff/eset_Mad2PB.RData'),
    opt('h', 'highlight', 'yaml', 'sets.yaml'),
    opt('o', 'outfile', 'RData', 'sets_mad2pb.RData'),
    opt('p', 'plotfile', 'pdf', 'sets_mad2pb.pdf'),
    arg('genesets', 'set expr RData files', arity='*',
        list.files("gsva_mad2pb", "\\.RData$", full.names=TRUE))
)

dset = io$load(args$dset)
highlight = io$read_yaml(args$highlight)

switch(args$highlight,
    "sets.yaml" = {
        use = highlight$expr_sets
        sets = io$load(args$genesets)[names(use)]
        expr = lapply(names(sets), function(s) sets[[s]][use[[s]],,drop=FALSE]) %>%
            narray::stack(along=1)
    },
    "genes.yaml" = {
        if (grepl("MILE", args$dset)) {
            highlight$genes = idmap$orthologue(highlight$genes, from="external_gene_name",
                    to="hgnc_symbol", dset="mmusculus_gene_ensembl") %>%
                unname() %>% na.omit()
            expr = dset$expr[intersect(highlight$genes, rownames(dset$expr)),]
        } else
            expr = dset$vs[highlight$genes,]
    },
    {
        stop("need to add highlight handler for ", args$highlight)
    }
)

switch(basename(args$dset),
    "eset_Mad2PB.RData" = {
        meta = as.data.frame(SummarizedExperiment::colData(dset$eset))
        groups = as.matrix(meta[c("Tcell", "Myeloid", "Other")])
    },
    "eset_MILE.RData" = {
        # copied from ../diff_expr/de_MILE.r
        keep = !is.na(dset$meta$type)
        expr = expr[,keep]
        groups = cbind(narray::mask(dset$meta$type[keep]),
            Hyperdip = (dset$meta$annot[keep] == "ALL with hyperdiploid karyotype")) + 0
        groups[,"Hyperdip"][groups[,"B_like"] == 0] = NA
    #        aneuploidy = pmin(dset$meta$aneuploidy[keep], 0.25)
    },
    {
        stop("need to add data set handler for ", args$dset)
    }
)

stopifnot(rownames(groups) == colnames(expr))

tmat = sapply(colnames(groups), function(g) {
    expr[,!is.na(groups[,g]) & groups[,g] == 1]
}, simplify=FALSE)

ecmp = lapply(tmat, reshape2::melt) %>%
    dplyr::bind_rows(.id="type") %>%
    dplyr::rename(set=Var1, sample=Var2, expr=value) %>%
    group_by(set) %>%
    mutate(z_expr = scale(expr)) %>%
    ungroup()

pdf(args$plotfile, 20, 15)
ggplot(ecmp, aes(x=z_expr, y=set)) +
    ggridges::geom_density_ridges(aes(fill=type), alpha=0.95)
dev.off()

save(expr, groups, file=args$outfile)

#if (grepl("rnaseq/assemble.RData", args$expr)) {
#    types = dset$idx$type
#    expr = dset$expr[match(c("Ets1", "Erg"), dset$genes),]
#    rownames(expr) = c("Ets1", "Erg")
#    meta = io$load(args$meta)
#    aneup0.3 = pmin(meta$aneuploidy[match(colnames(expr), meta$sample)], 0.3)
#    expr = rbind(expr, aneup0.3)
#} else if (grepl("E-GEOD", args$expr)) {
#    types = Biobase::pData(dset)$FactorValue..LEUKEMIA.CLASS.
#    expr = Biobase::exprs(dset)[c("ENSG00000134954", "ENSG00000157554"),]
#    rownames(expr) = c("ETS1", "ERG")
#} else {
#}
#
#mat = narray::stack(c(sets, list(expr)), along=1)[,colnames(expr)]
#if (is.character(types)) {
#    tmat = narray::split(mat, along=2, subsets=types)
#    mat2 = rbind(mat, narray::mask(types, along=1) + 0)
#} else { # mile2
#    tmat = sapply(colnames(types), function(t) {
#        mat[,!is.na(types[,t]) & types[,t] == 1]
#    }, simplify=FALSE)
#    mat2 = rbind(mat, t(types[,1:3]))
#}
