library(dplyr)
library(ggplot2)
sys = import('sys')
util = import('./util')

#' Fit tissue of origin (vs. mean) and pan-aneuploidy (tissue and type-corrected)
de_type = function(eset, ats=c("Tcell", "Myeloid", "Other")) {
    res = util$do_wald(eset, ~ tissue + type + aneuploidy, ex="tissue|aneuploidy")
    cancer_type = function(term) {
        DESeq2::design(eset) = formula(paste("~ tissue + ", term))
        DESeq2::estimateDispersions(eset) %>%
            DESeq2::nbinomWaldTest(maxit=1000) %>%
            util$extract_coef(term)
    }
    c(res, sapply(ats, cancer_type, simplify=FALSE))
}

#' Fit aneuploidy within cancer types
de_per_type = function(eset, ats=c("Tcell", "Myeloid", "Other")) {
    aneup_tissue = function(type) {
        eset = eset[,eset[[type]] == 1]
        if (length(unique(eset$tissue)) == 1)
            DESeq2::design(eset) = formula(paste("~ aneuploidy"))
        else
            DESeq2::design(eset) = formula(paste("~ aneuploidy + tissue"))
        DESeq2::estimateDispersions(eset) %>%
            DESeq2::nbinomWaldTest(maxit=1000) %>%
            util$extract_coef("aneuploidy")
    }
    setNames(lapply(ats, aneup_tissue), paste0(ats, ":aneuploidy"))
}

sys$run({
    args = sys$cmd$parse(
        opt('e', 'eset', 'gene expression rds', 'eset_Mad2PB.rds'),
        opt('c', 'config', 'yaml', '../config.yaml'),
        opt('n', 'network', 'rds', '../data/networks/E-GEOD-13159.rds'), # ignored
        opt('o', 'outfile', 'results rds', 'de_Mad2PB.rds'),
        opt('p', 'plotfile', 'pdf', 'de_Mad2PB.pdf')
    )

    hl = yaml::read_yaml(args$config)$highlight_de
    eset = readRDS(args$eset)$eset

    de_type = de_type(eset)
    de_per_type = de_per_type(eset)
    res = c(de_type, de_per_type)

    pdf(args$plotfile)
    for (rname in names(res))
        print(util$plot_volcano(res[[rname]], hl) + ggtitle(rname))
    dev.off()

    saveRDS(res, file=args$outfile)
})
