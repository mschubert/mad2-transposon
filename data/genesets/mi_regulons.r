library(dplyr)
b = import('base')
io = import('io')
sys = import('sys')
st = import('stats')
idmap = import('process/idmap')
aracne = import('tools/aracne')

args = sys$cmd$parse(
    opt('e', 'expr', 'expression RData', '../arrayexpress/E-GEOD-13159.RData'),
    opt('n', 'network', 'aracne results RData', '../networks/E-GEOD-13159.RData'),
    opt('o', 'outfile', 'gene set RData', 'MI_regulons.RData'))

expr = Biobase::exprs(io$load(args$expr))
rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
expr = expr[!is.na(rownames(expr)),]

get_cor = function(r, t) cor(expr[r,], expr[t,]) %catch% NA

regs = io$load(args$network) %>%
    mutate(cor = purrr::map2_dbl(Regulator, Target, get_cor),
           z = st$cor$fisher_r2z(cor),
           sign = ifelse(z > 0, "up", "down")) %>%
    filter(abs(z) > 0.5) %>%
    mutate(Regulator = paste(Regulator, sign, sep="_"))

regs = unstack(regs[c("Target", "Regulator")])

save(regs, file=args$outfile)
