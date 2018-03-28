library(dplyr)
b = import('base')
io = import('io')
st = import('stats')
idmap = import('process/idmap')
gs = import('data/genesets')
plt = import('plot')
df = import('data_frame')

# load samples, gene expression
samples = io$load('aneuploidy_mad2.RData') %>%
    mutate(tissue = b$grep("(S|T)", sample),
           mad2_class = factor(mad2_class, levels=c("low", "high")),
           aneup_class = factor(aneup_class, levels=c("low", "high")))
ploidy = io$load('compare_rna-scWGS_ploidy/ploidy_eT.RData')$ploidy
colnames(ploidy) = paste0(b$grep("([0-9]+)", colnames(ploidy)),
                          toupper(b$grep("(s|t|S|T)", colnames(ploidy))))
dset = io$load('../data/rnaseq/assemble.RData')
expr = dset$expr[rowSums(dset$counts) >= 10,]
colnames(expr) = paste0(b$grep("([0-9]+)", colnames(expr)),
                        toupper(b$grep("(s|t|S|T)", colnames(expr))))
narray::intersect(expr, ploidy, samples$sample, along=2)

# create gene and ploidy matrices
genes = t(idmap$gene(expr, to="external_gene_name", dset="mmusculus_gene_ensembl"))
gene_chrs = idmap$gene(colnames(genes), from="external_gene_name",
                       to="chromosome_name", dset="mmusculus_gene_ensembl")
valid = gene_chrs %in% as.character(1:19) #TODO: 'X' missing in ploidy obj
genes = genes[,names(gene_chrs)[valid]]
ploidy = t(ploidy[gene_chrs[valid],])

# go annotations from ensembl, scores matrix
go = gs$go(dset="mmusculus_gene_ensembl", genes="external_gene_name") %>%
    filter(external_gene_name %in% colnames(genes)) %>%
    transmute(gene = external_gene_name,
              set = paste(id, name)) %>%
    unstack()
size = sapply(go, length)
go = go[size >= 5 & size <= 200]
scores = GSVA::gsva(t(genes), go, kernel=FALSE, parallel.sz=1) %>% t()

# gene expression associations for mad2 and aneuploidy classes
assocs = function(mat, subs, include) {
    mat = mat[include,]
    blocking = samples$tissue[include]
    subs = subs[include]
    st$lm(mat ~ blocking + subs) %>%
        filter(!grepl("blocking", term)) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
}

assocs2 = function(mat, subs, include) {
    mat = mat[include,]
    ploidy = ploidy[include,]
    blocking = samples$tissue[include]
    subs = subs[include]

    lapply(1:ncol(mat), function(n) {
            lm(mat[,n] ~ blocking + ploidy[,n] + subs) %>% broom::tidy()
        }) %>%
        setNames(colnames(mat)) %>%
        df$bind_rows("mat") %>%
        filter(grepl("subs", term)) %>%
        select(-term) %>%
        mutate(size=sum(include),
               adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value)
}

# removing thymus here makes a tiny difference, basically the same genes
gene_mad2 = assocs(genes, samples$mad2_class, samples$aneup_class == "low")
set_mad2 = assocs(scores, samples$mad2_class, samples$aneup_class == "low")
gene_aneup = assocs(genes, samples$aneup_class, samples$mad2_class == "low")
set_aneup = assocs(scores, samples$aneup_class, samples$mad2_class == "low")

gene_mad2_cor = assocs2(genes, samples$mad2_class, samples$aneup_class == "low")
#set_mad2_cor = assocs2(scores, samples$mad2_class, samples$aneup_class == "low")
gene_aneup_cor = assocs2(genes, samples$aneup_class, samples$mad2_class == "low")
#set_aneup_cor = assocs2(scores, samples$aneup_class, samples$mad2_class == "low")

save(gene_mad2, set_mad2, gene_aneup, set_aneup,
     gene_mad2_cor, gene_aneup_cor, file="diff_expr_classes.RData")

# volcano plots
volcano = function(assocs, ...) {
    assocs %>%
        mutate(label = mat) %>%
        plt$p_effect("adj.p", thresh=0.1) %>%
        plt$volcano(..., repel=TRUE, label_top=50)
}

pdf("diff_expr_classes.pdf")
volcano(gene_mad2 %>% filter(mat != "Mad2l1")) + ggtitle("Mad2 high, genes")
volcano(gene_mad2_cor %>% filter(mat != "Mad2l1")) + ggtitle("Mad2 high, genes - corrected for ploidy")
volcano(set_mad2, text.size=1.5) + ggtitle("Mad2 high, gene sets")
volcano(gene_aneup) + ggtitle("Aneuploidy high, genes")
volcano(gene_aneup_cor) + ggtitle("Aneuploidy high, genes - corrected for ploidy")
volcano(set_aneup, text.size=1.5) + ggtitle("Aneuploidy high, gene sets")
dev.off()
