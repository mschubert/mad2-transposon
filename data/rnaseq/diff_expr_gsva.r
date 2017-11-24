library(DESeq2)
b = import('base')
io = import('io')
st = import('stats')
gs = import('data/genesets')
idmap = import('process/idmap')

if (is.null(module_name())) {
    # load the gene expression data + sample types
    dset = io$load('assemble.RData')
    expr = dset$expr[rowSums(dset$counts) >= 10,]
    index = dset$idx %>%
        transmute(id = id,
                  tissue = tissue,
                  type = `Tumour type`,
                  mad2 = `Mad2 levels`,
                  stage = `Early/ Late`)
    index$tissue[is.na(index$tissue)] = "other"

    # load Gene Ontology sets from Ensembl
    go = gs$go(dset="mmusculus_gene_ensembl", genes="ensembl_gene_id") %>%
        filter(ensembl_gene_id %in% rownames(expr)) %>%
        transmute(gene = ensembl_gene_id,
                  set = paste(id, name)) %>%
        unstack()
    size = sapply(go, length)
    go = go[size >= 5 & size <= 200]

    # no need to scale expr here, GSVA does it
    scores = GSVA::gsva(expr, go, kernel=FALSE, parallel.sz=1) %>% t()

    cond = cbind(intercept = rep(1, nrow(index)),
                 thymus_over_spleen = index$tissue == "thymus",
                 other_over_spleen = index$tissue == "other",
                 myeloid_leukemia = index$type == "Myeloid leukaemia",
                 stage_late = index$stage == "Late",
                 mad2_high = index$mad2 == "High") + 0
    rownames(cond) = index$id

    process = . %>%
        mutate(term = sub("^cond", "", term)) %>%
        filter(term != "intercept") %>%
        select(-size) %>%
        group_by(term) %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        ungroup() %>%
        arrange(adj.p, p.value) %>%
        filter(adj.p < 0.2)

    gene = t(expr)
    assocs_gene = st$lm(gene ~ 0 + cond, atomic="cond") %>% process() %>%
        mutate(gene = idmap$gene(gene, to="external_gene_name", dset="mmusculus_gene_ensembl"))
    assocs_sets = st$lm(scores ~ 0 + cond, atomic="cond") %>% process()

    write.table(assocs_gene, file="de_genes.tsv", row.names=FALSE, quote="", sep="\t")
    write.table(assocs_sets, file="de_sets.tsv", row.names=FALSE, quote="", sep="\t")
}
