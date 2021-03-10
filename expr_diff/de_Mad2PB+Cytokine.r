library(dplyr)
library(ggplot2)
library(DESeq2)
sys = import('sys')
plt = import('plot')
gset = import('data/genesets')
util = import('./util')

args = sys$cmd$parse(
    opt('e', 'eset', 'gene expression rds', 'eset_Mad2PB+Cytokine.rds'),
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('n', 'network', 'rds', '../data/networks/E-GEOD-13159.rds'), # ignored
    opt('o', 'outfile', 'results rds', 'de_Mad2PB+Cytokine.rds'),
    opt('p', 'plotfile', 'pdf', 'de_Mad2PB+Cytokine.pdf'),
    arg('sets', 'gene set .rds', arity='*',
        list.files("../data/genesets/mouse", "\\.rds", full.names=TRUE))
)

hl = c('Ccl5', 'Cxcl10', 'Ifna1', 'Ifnb1', 'Tnf', 'Il6', 'Ccl20', 'Cxcl1', 'Ccl2', 'Il2', 'Ifng',
    'Ccr2', 'Ccr5', 'Cxcr3', 'Ifnar1', 'Ifnar2', 'Tlr3', 'Tlr4', 'Tlr7', 'Tlr8', 'Tlr9',
    'Il10ra', 'Il6ra', 'Ccr1', 'Ccr4', 'Tnfrsf1a', 'Tradd', 'Cxcr2', 'Ifngr1', 'Ifngr2',
    'Samd9l', 'Isg15', 'Ifi27l2a', 'Ifit1', 'Oas1a', 'Oas1b', 'Isg54', 'Isg56')

gset = readRDS(grep("HALLMARK", args$sets, value=TRUE))
gset = unlist(gset[c('HALLMARK_INTERFERON_GAMMA_RESPONSE',
              'HALLMARK_INFLAMMATORY_RESPONSE', 'HALLMARK_TNFA_SIGNALING_VIA_NFKB')])
gset = c(gset, hl)
eset = readRDS(args$eset)$eset
eset = eset[intersect(gset, rownames(eset)),]

warning("Not in eset: ", paste(hl[!hl %in% rownames(eset)], collapse=", "))

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

pdf(args$plotfile)
for (rname in names(res)) {
    message(rname)
    print(util$plot_volcano(res[[rname]], hl) + ggtitle(rname))
}
dev.off()

saveRDS(res, file=args$outfile)
