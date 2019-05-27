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
    opt('e', 'eset', 'gene expression RData', 'eset_Mad2PB.RData'),
    opt('c', 'config', 'yaml', '../config.yaml'),
    opt('n', 'network', 'RData', '../data/networks/E-GEOD-13159.RData'), # ignored
    opt('o', 'outfile', 'results RData', 'de_Mad2PB.RData'),
    opt('p', 'plotfile', 'pdf', 'de_Mad2PB.pdf'),
    arg('sets', 'gene set .RData', arity='*',
        list.files("../data/genesets/mouse", "\\.RData", full.names=TRUE)))

eset = io$load(args$eset)
ee = eset$vs
eset = eset$eset

# add stat1_act as col to eset (TODO: remove genes used for activity?
stat1 = ee["Stat1",]
rest = ee[-which(rownames(ee) == "Stat1"),]
cc = data.frame(
    gene = rownames(rest),
    cor = narray::map(rest, along=2, function(x) cor(x, stat1))
) %>% arrange(-cor)
stat1act = GSVA::gsva(ee, list(Stat1_act = head(cc$gene, 500)))
df = colData(eset)
df$stat1act = stat1act[1,rownames(df)]
colData(eset) = df

# fit tissue of origin and pan-aneuploidy
res = util$do_wald(eset, ~ tissue + type + aneuploidy * stat1act, ex="aneuploidy|stat1act")
## fit cancer type vs mean
#cancer_type = function(term) {
#    design(eset) = formula(paste("~ tissue + ", term))
#    DESeq2::estimateDispersions(eset) %>%
#        DESeq2::nbinomWaldTest(maxit=1000) %>%
#        util$extract_coef(term)
#}
#ats = c("Tcell", "Myeloid", "Other")
#res = c(res, sapply(ats, cancer_type, simplify=FALSE))

## fit aneuploidy within cancer types
#aneup_tissue = function(type) {
#    eset = eset[,eset[[type]] == 1]
#    if (length(unique(eset$tissue)) == 1)
#        design(eset) = formula(paste("~ aneuploidy * stat1act"))
#    else
#        design(eset) = formula(paste("~ tissue + aneuploidy * stat1act"))
#    DESeq2::estimateDispersions(eset) %>%
#        DESeq2::nbinomWaldTest(maxit=1000) %>%
#        util$extract_coef("stat1act")
#}
#res = c(res, setNames(lapply(ats, aneup_tissue), paste0(ats, ":stat1act")))

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

save(res, file=args$outfile)
