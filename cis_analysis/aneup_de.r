io = import('io')
sys = import('sys')
idmap = import('process/idmap')
plt = import('plot')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression', '../data/rnaseq/assemble.RData'),
    opt('c', 'cis', 'CIS assocs', 'poisson.RData'),
    opt('a', 'aneup', 'aneuploidy RData', '../ploidy_compare/analysis_set.RData'),
    opt('o', 'outfile', '.RData', 'aneup_de.RData'),
    opt('p', 'plotfile', 'pdf', 'aneup_de.pdf'))

genes = io$load(args$cis)$result %>%
    filter(adj.p < 1e-3)

aneup = io$load(args$aneup)

dset = io$load(args$expr)
expr = dset$expr %>%
    idmap$gene(to="external_gene_name", dset="mmusculus_gene_ensembl")
expr = expr[genes$external_gene_name,]

narray::intersect(expr, aneup$sample, along=2)
aneup = aneup$aneup
models = narray::lambda(~ lm(expr ~ aneup), along=c(expr=1))
result = models %>%
    mutate(ins_p = genes$p.value[match(expr, genes$external_gene_name)],
           result = lapply(result, broom::tidy)) %>%
    tidyr::unnest() %>%
    filter(term == "aneup") %>%
    select(-term) %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(adj.p, p.value)

pdf(args$plotfile)
result %>%
    mutate(label = expr,
           size = -log2(ins_p)) %>%
    plt$color$p_effect("p.value") %>%
    plt$volcano()

for (g in result %>% filter(adj.p < 0.1) %>% pull(expr)) {
    i = which(g == models$expr)
    plot(model.frame(models$result[[i]])[,2:1], main=models$expr[i])
}
dev.off()

save(result, file=args$outfile)
