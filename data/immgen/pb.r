library(dplyr)
io = import('io')
sys = import('sys')
dset = import('./dset')

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'dset.rds'),
    opt('o', 'outfile', 'rds', 'pb.rds'),
    opt('p', 'plotfile', 'pdf', 'pb.pdf')
)

dset = readRDS(args$infile)
meta = dset$meta %>%
    filter(!is.na(annot))
expr = dset$expr
narray::intersect(expr, meta$fname, along=2)

## select informative genes for lineages
gene_variances = matrixStats::rowVars(expr)
hivar = head(order(gene_variances, decreasing=TRUE), 1000)
expr = expr[hivar,]
#x = readr::read_tsv("mart_export.txt") %>%
#    filter(`GO term accession` == "GO:0003700") %>%
#    pull(`Gene stable ID`) %>%
#    intersect(rownames(expr))
#expr = expr[x,]

# batch-correct mouse tumors
tps = readRDS("../rnaseq/assemble.rds")
tpse = tps$expr
tpsm = tps$idx

e2 = narray::stack(expr, tpse, along=2)
e2 = e2[rowSums(is.na(e2)) == 0,]
merge = sva::ComBat(e2, batch=c(rep("immgen", ncol(expr)), rep("tps", ncol(tpse))), par.prior=TRUE)
expr = merge[,1:ncol(expr)]
tpse = merge[,-c(1:ncol(expr))]

# sanity checks for batch correction
#dset$plot_pca(meta, expr, color)

# train LDA model on reference set
library(MASS)
mod = lda(meta$annot ~ ., data=as.data.frame(t(expr)))


# predict label for each mouse tumor
res = predict(mod, newdata=as.data.frame(t(tpse)))
tpsm %>% cbind(pred=res$class) %>% dplyr::select(mouse, genotype, tissue, type, pred) %>% as.data.frame()


# plot mouse tumors along differentiation tree
