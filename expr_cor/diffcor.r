library(dplyr)
library(corrplot)
library(cowplot)
library(ggraph)
io = import('io')
sys = import('sys')
st = import('stats')
plt = import('plot')

cor_diff_2d = function(mlist) {
    nrm = function(x) do.call(narray::melt, x)
    dt = function(x, y) nrm(st$cor$diff_test(t(mlist[[x]]), t(mlist[[y]])))
    dct = function(x) nrm(list(cor=cor(t(mlist[[x]])), p.value=st$cor$test(t(mlist[[x]]))))

    cors = data.frame(type1=names(mlist), type2=names(mlist), stringsAsFactors=FALSE) %>%
        mutate(res = purrr::map(type1, dct)) %>%
        tidyr::unnest()
        
    diffs = expand.grid(names(mlist), names(mlist), stringsAsFactors=FALSE) %>%
        dplyr::rename(type1 = Var1, type2 = Var2) %>%
        filter(type1 < type2) %>%
        mutate(res = purrr::map2(type1, type2, dt)) %>%
        tidyr::unnest() %>%
        dplyr::rename(cor = delta_cor)

    both = dplyr::bind_rows(cors, diffs) %>%
        mutate(label = -log10(p.value),
               label = ifelse(is.infinite(label), 0, label)) %>%
        group_by(type1, type2) %>%
        mutate(rank = rank(-label),
               label = ifelse(rank <= 10, label, NA)) %>%
        ungroup() %>%
        mutate(label = floor(label))
    ord = both %>%
        group_by(Var1, Var2) %>%
        summarize(sumcor = sum(cor)) %>%
        plt$cluster(sumcor ~ Var1 + Var2)
    both$Var1 = factor(both$Var1, levels=levels(ord$Var1))
    both$Var2 = factor(both$Var2, levels=levels(ord$Var2))
    ggplot(both, aes(x=Var1, y=Var2)) +
        facet_grid(type1 ~ type2) +
        geom_tile(aes(size=abs(cor), fill=cor)) + #, alpha=abs(cor))) +
#        ggforce::geom_circle(aes(x0=Var1, y0=Var2, r=1)) +
        scale_fill_gradient2(low="red", mid="white", high="blue") +
        coord_fixed() +
        geom_text(aes(label=label)) +
        theme(axis.text.x = element_text(angle=90, hjust=1))
}

args = sys$cmd$parse(
    opt('e', 'expr', 'expr RData', '../expr_sets/genes_mile.RData'),
    opt('p', 'plotfile', 'pdf', 'diffcor-genes_mile.pdf')
)

dset = io$load(args$expr)
mat = dset$expr
tmat = sapply(colnames(dset$groups), function(t) {
    dset$expr[,!is.na(dset$groups[,t]) & dset$groups[,t] == 1]
}, simplify=FALSE)

pdf(16, 14, file=args$plotfile)
cor_diff_2d(tmat)
dev.off()
