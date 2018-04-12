library(DESeq2)
b = import('base')
io = import('io')
st = import('stats')
plt = import('plot')
gs = import('data/genesets')
idmap = import('process/idmap')

assocs_genes = function(expr) {
    # load Gene Ontology sets from Ensembl
    gene = t(expr)
    assocs_gene = st$lm(gene ~ type + aneuploidy, data=index) %>%
        mutate(gene = idmap$gene(gene, to="external_gene_name", dset="mmusculus_gene_ensembl")) %>%
        arrange(p.value)
}

assocs_sets = function(expr) {
    gene = t(expr)

    go = gs$go(dset="mmusculus_gene_ensembl", genes="ensembl_gene_id") %>%
        filter(ensembl_gene_id %in% rownames(expr)) %>%
        transmute(gene = ensembl_gene_id,
                  set = paste(id, name)) %>%
        unstack()
    size = sapply(go, length)
    go = go[size >= 5 & size <= 200]

    # no need to scale expr here, GSVA does it
    scores = GSVA::gsva(expr, go, kernel=FALSE, parallel.sz=1) %>% t()
    assocs_sets = st$lm(scores ~ type + aneuploidy, data=index) %>%
        arrange(p.value)
}

volcano = function(assocs, cond, repel=TRUE, label_top=50, ...) {
    assocs %>%
        filter(term == cond) %>%
        mutate(adj.p=p.adjust(p.value, method="fdr"), label = scores, size=5) %>%
        plt$p_effect("adj.p", thresh=0.1) %>%
        plt$volcano(..., repel=repel, label_top=label_top)
}

if (is.null(module_name())) {
    # load the gene expression data + sample types
    aneuploidy = io$load('../aneuploidy/ploidy_from_rnaseq.RData')$aneuploidy
    dset = io$load('../data/rnaseq/assemble.RData')
    expr = dset$expr[rowSums(dset$counts) >= 10,]
    index = dset$idx %>%
        transmute(id = id,
                  tissue = tissue,
                  type = `Tumour type`,
                  mad2 = `Mad2 levels`,
                  stage = `Early/ Late`) %>%
        full_join(data.frame(id=names(aneuploidy), aneuploidy=unname(aneuploidy)), by = "id")
    index$tissue[is.na(index$tissue)] = "other" # should this be spleen?

    #TODO?: abs(scale(expr)) if parts of set up+down
    assocs_gene = assocs_genes(expr)
    assocs_sets = assocs_sets(expr)

    write.table(assocs_gene, file="de_genes.tsv", row.names=FALSE, quote=FALSE, sep="\t")
    write.table(assocs_sets, file="de_sets.tsv", row.names=FALSE, quote=FALSE, sep="\t")

    assocs_gene$scores = assocs_gene$gene
    pdf("volcano_genes.pdf", width=10, height=10)
    for (subs in unique(assocs_gene$term))
        print(volcano(assocs_gene, subs) + ggtitle(subs))
    dev.off()

    pdf("volcano_sets.pdf", width=10, height=10)
    for (subs in unique(assocs_sets$term))
        print(volcano(assocs_sets, subs, label_top=20) + ggtitle(subs))
    dev.off()

    #TODO: network enrichment plot(?)
}
