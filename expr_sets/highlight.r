library(dplyr)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('d', 'dset', 'expr RData', '../expr_diff/eset_Mad2PB.RData'),
    opt('h', 'highlight', 'yaml', 'interesting.yaml'),
    opt('o', 'outfile', 'RData', 'hl_mad2pb.RData'),
    arg('genesets', 'RData files', arity='*',
        list.files("../expr_sets/gsva_mad2pb", "\\.RData$", full.names=TRUE))
)

dset = io$load(args$dset)
highlight = io$read_yaml(args$highlight)$expr_sets
sets = io$load(args$genesets)
sets = lapply(names(highlight), function(s) sets[[s]][highlight[[s]],,drop=FALSE])

if (grepl("rnaseq/assemble.RData", args$expr)) {
    types = dset$idx$type
    expr = dset$expr[match(c("Ets1", "Erg"), dset$genes),]
    rownames(expr) = c("Ets1", "Erg")
    meta = io$load(args$meta)
    aneup0.3 = pmin(meta$aneuploidy[match(colnames(expr), meta$sample)], 0.3)
    expr = rbind(expr, aneup0.3)
} else if (grepl("E-GEOD", args$expr)) {
    types = Biobase::pData(dset)$FactorValue..LEUKEMIA.CLASS.
    expr = Biobase::exprs(dset)[c("ENSG00000134954", "ENSG00000157554"),]
    rownames(expr) = c("ETS1", "ERG")
} else {
    # copied from ../diff_expr/de_MILE.r
    keep = !is.na(dset$meta$type)
    expr = dset$expr[c("ETS1", "ERG"),keep]
    types = cbind(narray::mask(dset$meta$type[keep]),
        Hyperdip = (dset$meta$annot[keep] == "ALL with hyperdiploid karyotype")) + 0
    types[,"Hyperdip"][types[,"B_like"] == 0] = NA
#        aneuploidy = pmin(dset$meta$aneuploidy[keep], 0.25)
}

mat = narray::stack(c(sets, list(expr)), along=1)[,colnames(expr)]
if (is.character(types)) {
    tmat = narray::split(mat, along=2, subsets=types)
    mat2 = rbind(mat, narray::mask(types, along=1) + 0)
} else { # mile2
    tmat = sapply(colnames(types), function(t) {
        mat[,!is.na(types[,t]) & types[,t] == 1]
    }, simplify=FALSE)
    mat2 = rbind(mat, t(types[,1:3]))
}

save(sets, groups, file=args$outfile)
