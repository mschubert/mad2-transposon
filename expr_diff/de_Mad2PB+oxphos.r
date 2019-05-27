library(dplyr)
library(ggplot2)
library(DESeq2)
io = import('io')
sys = import('sys')
plt = import('plot')
idmap = import('process/idmap')
gset = import('data/genesets')
util = import('./util')

args = sys$cmd$parse(
    opt('e', 'eset', 'gene expression RData', 'eset_Mad2PB+MycV1.RData'),
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('n', 'network', 'RData', '../data/networks/E-GEOD-13159.RData'), # ignored
    opt('o', 'outfile', 'results RData', 'de_Mad2PB+MycV1.RData'),
    opt('p', 'plotfile', 'pdf', 'de_Mad2PB+MycV1.pdf'),
    arg('sets', 'gene set .RData', arity='*',
        list.files("../data/genesets/mouse", "\\.RData", full.names=TRUE)))

gset = io$load(grep("HALLMARK", args$sets, value=TRUE))$HALLMARK_OXIDATIVE_PHOSPHORYLATION
eset = io$load(args$eset)$eset
eset = eset[intersect(gset, rownames(eset)),]

# fit tissue of origin and pan-aneuploidy
res = util$do_wald(eset, ~ tissue + type + aneuploidy, ex="tissue|aneuploidy")
# fit cancer type vs mean
cancer_type = function(term) {
    design(eset) = formula(paste("~ tissue + ", term))
    DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomWaldTest(maxit=1000) %>%
        util$extract_coef(term)
}
ats = c("Tcell", "Myeloid", "Other")
res = c(res, sapply(ats, cancer_type, simplify=FALSE))

# fit aneuploidy within cancer types
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
res = c(res, setNames(lapply(ats, aneup_tissue), paste0(ats, ":aneuploidy")))

hl = io$read_yaml(args$config)$highlight_de

pdf(args$plotfile)
for (rname in names(res)) {
    message(rname)
    print(util$plot_volcano(res[[rname]], hl) + ggtitle(rname))
}
dev.off()

save(res, file=args$outfile)
