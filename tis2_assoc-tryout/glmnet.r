# beta regression
#
# use the fraction of reads that are in transposon
# aneuploidy scores as continuous variable
# test whether or not each gene has insertion w/ increasing aneuploidy score
library(dplyr)
library(glmnet)
io = import('io')
plt = import('./plot')
sys = import('sys')

do_fit = function(data) {
    x = narray::construct(reads ~ sample + gene_name, data=data, fill=0)
#    x = x[,narray::map(x, along=1, sd) > 0]
    y = data %>% select(sample, aneup) %>% distinct()
    narray::intersect(x, y$sample, along=1)
    ll = list(cv = cv.glmnet(x, y$aneup, alpha=1))
    ll$opt = glmnet(x, y$aneup)
    ll$feats = as.matrix(coef(ll$opt, s=ll$cv$lambda.min))[,1]
    ll$feats = ll$feats[ll$feats != 0 & names(ll$feats) != "(Intercept)"]
    ll
}

args = sys$cmd$parse(
    opt('s', 'samples', 'min samples w/ insertions in gene', '5'),
    opt('d', 'decay', 'decay factor for flanking genes', '0.7'),
    opt('o', 'outfile', 'file to save to', 'glmnet.RData'),
    opt('p', 'plotfile', 'pdf', 'glmnet.pdf'))

cis = io$load("dset.RData") %>%
    filter(!is.na(aneup))

dset = . %>%
    group_by(sample) %>%
    mutate(reads = log(reads),
           reads = as.numeric(args$decay)^hit_dist * reads / sum(reads[hit_dist == 0])) %>%
    group_by(sample, gene_name, aneup) %>%
    summarize(reads = sum(reads)) %>%
    ungroup() %>%
    group_by(gene_name) %>%
    filter(n() >= as.integer(args$samples)) %>%
    ungroup()

subs = list(
    hits = dset(filter(cis, hit_dist == 0)),
    hits_cancer = dset(filter(cis, hit_dist == 0, known_cancer)),
    near = dset(cis),
    near_cancer = dset(filter(cis, known_cancer))
)
fits = lapply(subs, do_fit)

pdf(args$plotfile)
for (i in seq_along(subs)) {
    plot(fits[[i]]$cv)
    plot(fits[[i]]$opt)
    print(plt$plot_tiles(subs[[i]], genes=names(fits[[i]]$feats)))
}
dev.off()

save(fits, subs, file=args$outfile)
