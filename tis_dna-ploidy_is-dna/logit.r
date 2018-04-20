# logistic regression
#
# 20(ish) read cutoff per transposon insertion
# aneuploidy scores as continuous variable
# test whether or not each gene has insertion w/ increasing aneuploidy score
library(dplyr)
io = import('io')

dset = io$load("dset.RData")

do_fit = function(i, obj="cis_hits", min_reads=20) {
    gene = dset[[obj]][i,] >= min_reads
    aneup = dset$aneup$aneup
    glm(gene ~ aneup, family=binomial(link='logit')) %>%
        broom::tidy() %>%
        filter(term == "aneup") %>%
        select(-term)
}
result = sapply(rownames(dset$cis_hits), do_fit, simplify=FALSE) %>%
    dplyr::bind_rows(.id="gene") %>%
    arrange(p.value)
