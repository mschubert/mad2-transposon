library(dplyr)
io = import('io')
sys = import('sys')
gset = import('data/genesets')
util = import('../expr_diff/util')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('s', 'sigfile', 'rds', 'sig_statKO.rds'),
    opt('e', 'expr', 'RData', '../expr_diff/eset_Mad2PB.RData'),
    opt('p', 'plotfile', 'pdf', 'de_statKO.pdf'),
    arg('sets', 'RDatas with gene sets', arity='*',
        list.files("../data/genesets/mouse", full.names=TRUE)))

sig = readRDS(args$sigfile)
dset = io$load(args$expr)
eset = dset$eset

# add signature score to eset colData
if (is.numeric(sig)) { # z-score vector
    vs = dset$vs
    narray::intersect(sig, vs, along=1)
    scores = (t(as.matrix(sig)) %*% vs)[1,]
} else { # gene symbol set
    scores = GSVA::gsva(dset$vs, list(sig=sig))[1,]
}
colData(eset)$scores = scores

# fit tissue of origin and pan-aneuploidy
res = list(aneuploidy=util$do_wald(eset, ~ tissue + type + scores + aneuploidy, ex="aneuploidy"))

# fit aneuploidy within cancer types TODO: do we even care about this?
aneup_tissue = function(type) {
    eset = eset[,eset[[type]] == 1]
    if (length(unique(eset$tissue)) == 1)
        design(eset) = formula(paste("~ aneuploidy"))
    else
        design(eset) = formula(paste("~ aneuploidy + tissue"))
    DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomWaldTest(maxit=1000) %>%
        util$extract_coef("aneuploidy")
}
#res = c(res, setNames(lapply(ats, aneup_tissue), paste0(ats, ":aneuploidy")))

sets = io$load(args$sets) %>%
    setNames(tools::file_path_sans_ext(basename(args$sets))) %>%
    lapply(function(x) gset$filter(x, min=5, valid=rownames(eset)))

hl = io$read_yaml(args$config)$highlight_de

pdf(args$plotfile)
for (rname in names(res)) {
    message(rname)
    print(util$plot_volcano(res[[rname]], hl) + ggtitle(rname))
    for (sname in names(sets)) {
        title = paste(rname, sname)
        message(title)
        print(util$plot_gset(res[[rname]], sets[[sname]]) + ggtitle(title))
    }
}
dev.off()
