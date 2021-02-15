library(dplyr)
io = import('io')
sys = import('sys')

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

# batch-correct mouse tumors
tps = io$load("../rnaseq/assemble.RData")
tpse = tps$expr
tpsm = tps$idx

e2 = narray::stack(expr, tpse, along=2)
e2 = e2[rowSums(is.na(e2)) == 0,]
merge = sva::ComBat(e2, batch=c(rep("immgen", ncol(expr)), rep("tps", ncol(tpse))), par.prior=TRUE)
expr = merge[,1:ncol(expr)]
tpse = merge[,-c(1:ncol(expr))]

# train LDA model on reference set
library(MASS)
mod = lda(meta$annot ~ ., data=as.data.frame(t(expr)))


# predict label for each mouse tumor

# plot mouse tumors along differentiation tree
