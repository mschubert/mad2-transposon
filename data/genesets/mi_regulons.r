library(dplyr)
b = import('base')
sys = import('sys')
st = import('stats')
idmap = import('process/idmap')
aracne = import('tools/aracne')

args = sys$cmd$parse(
    opt('e', 'expr', 'expression RData', '../arrayexpress/E-GEOD-13159.rds'),
    opt('n', 'network', 'aracne results RData', '../networks/E-GEOD-13159.rds'),
    opt('h', 'human', 'save to RData', 'human/KEA_2015.rds'),
    opt('m', 'mouse', 'save to RData', 'mouse/KEA_2015.rds')
)

expr = Biobase::exprs(readRDS(args$expr))
rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
expr = expr[!is.na(rownames(expr)),]

get_cor = function(r, t) cor(expr[r,], expr[t,]) %catch% NA
to_mouse = function(x)
    unname(idmap$orthologue(x, from="external_gene_name", to="mgi_symbol"))

regs = readRDS(args$network) %>%
    mutate(cor = purrr::map2_dbl(Regulator, Target, get_cor),
           z = st$cor$fisher_r2z(cor),
           sign = ifelse(z > 0, "up", "down"),
           Target_mouse = to_mouse(Target)) %>%
    filter(abs(z) > 0.5,
           !is.na(Regulator),
           !is.na(Target)) %>%
    mutate(Regulator = paste(Regulator, sign, sep="_"))

human = unstack(regs[c("Target", "Regulator")])
mouse = unstack(regs[c("Target_mouse", "Regulator")])

saveRDS(human, file=args$human)
saveRDS(mouse, file=args$mouse)
