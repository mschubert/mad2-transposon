library(dplyr)
b = import('base')
io = import('io')
st = import('stats')
idmap = import('process/idmap')
gs = import('data/genesets')
plt = import('plot')

# load samples, gene expression
samples = io$load('aneuploidy_mad2.RData') %>%
    mutate(tissue = b$grep("(S|T)", sample),
           mad2_class = factor(mad2_class, levels=c("low", "high")),
           aneup_class = factor(aneup_class, levels=c("low", "high")))
dset = io$load('../data/rnaseq/assemble.RData')
expr = dset$expr[rowSums(dset$counts) >= 10,]
colnames(expr) = paste0(b$grep("([0-9]+)", colnames(expr)),
                        toupper(b$grep("(s|t|S|T)", colnames(expr))))
narray::intersect(expr, samples$sample, along=2)

# create gene matrix
genes = t(idmap$gene(expr, to="external_gene_name", dset="mmusculus_gene_ensembl"))

# go annotations from ensembl, scores matrix
go = gs$go(dset="mmusculus_gene_ensembl", genes="ensembl_gene_id") %>%
    filter(ensembl_gene_id %in% rownames(expr)) %>%
    transmute(gene = ensembl_gene_id,
              set = paste(id, name)) %>%
    unstack()
size = sapply(go, length)
go = go[size >= 5 & size <= 200]
scores = GSVA::gsva(expr, go, kernel=FALSE, parallel.sz=1) %>% t()

# gene expression associations for mad2 and aneuploidy classes
assocs = function(mat, blocking, subs, include) {
    mat = mat[include,]
    blocking = blocking[include]
    subs = subs[include]
    st$lm(mat ~ blocking + subs) %>%
        filter(!grepl("blocking", term)) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
}

# removing thymus here makes a tiny difference, basically the same genes
gene_mad2 = assocs(genes, samples$tissue, samples$mad2_class, samples$aneup_class == "low")
set_mad2 = assocs(scores, samples$tissue, samples$mad2_class, samples$aneup_class == "low")
gene_aneup = assocs(genes, samples$tissue, samples$aneup_class, samples$mad2_class == "low")
set_aneup = assocs(scores, samples$tissue, samples$aneup_class, samples$mad2_class == "low")
save(gene_mad2, set_mad2, gene_aneup, set_aneup, file="diff_expr_classes.RData")

# volcano plots
volcano = function(assocs, ...) {
    assocs %>%
        mutate(label = mat) %>%
        plt$p_effect("adj.p", thresh=0.1) %>%
        plt$volcano(..., repel=TRUE, label_top=50)
}

pdf("diff_expr_classes.pdf")
volcano(gene_mad2 %>% filter(mat != "Mad2l1")) + ggtitle("Mad2 high")
volcano(set_mad2, text.size=1.5) + ggtitle("Mad2 high")
volcano(gene_aneup %>% filter(mat != "Mad2l1")) + ggtitle("Aneuploidy high")
volcano(set_aneup, text.size=1.5) + ggtitle("Aneuploidy high")
dev.off()
