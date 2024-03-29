library(DESeq2)
b = import('base')
io = import('io')
st = import('stats')
gs = import('data/genesets')

#' Compute differential expression of genes
#'
#' Can do e.g. DESeq2::plotMA() on result
#'
#' @param counts  A count matrix [genes x samples]
#' @param index   A data.frame with sample metadata
#' @param fml     A formula specifying the conditions
#'                result will be last term, rest regressed out
#' @return        A data.frame with fields: baseMean, log2FoldChange,
#'                lfcSE (shrunk log2 FC), pvalue, etc.
de_genes = function(counts, fml) {
    cds = DESeqDataSetFromMatrix(counts, index, fml)
    dds = DESeq2::DESeq(cds)
    res = DESeq2::results(dds)
#    DESeq2::lfcShrink(res, coef=2)
}

#' Compute differential expression of gene sets
#'
#' @param de_genes  The DE gene result data.frame from de_genes(...)
#' @param sets      A list of character vectors representing gene names
#' @return
de_sets = function(de_genes, sets) {
    deg = na.omit(de_genes[c('log2FoldChange', 'pvalue')])
    piano::runGSA(#geneLevelStats = setNames(de_genes$stat, rownames(de_genes)),
                  geneLevelStats = setNames(de_genes$pvalue, rownames(de_genes)),
                  directions = setNames(de_genes$stat, rownames(de_genes)),
#                  geneSetStat = "gsea",
#                  signifMethod = "geneSampling", # samplePermutation needs permStats obj
                  gsc = piano::loadGSC(stack(sets)))
}

#' Produce a tidy data.frame from the paino run
#'
#' @param res   Result object from piano
#' @param type  Which type of statistic to return; multiple allowed
#' @return      A data.frame
clean_piano = function(res, type=c("p", "stat", "adj.p")) {
    type = sub("p", "p[^A]", type)
    type = sub("adj.p", "pAdj", type)
    type = paste(paste0("^", type), collapse="|")
    fields = grep(type, names(gsa), value=TRUE)

    tibble::as_tibble(cbind(set=names(gsa$gsc), as.data.frame(gsa[fields])))
}

if (is.null(module_name())) {
    # load the gene expression data + sample types
    dset = io$load('assemble.RData')
    counts = dset$counts[rowSums(dset$counts) >= 10,]
    index = dset$idx %>%
        transmute(id = id,
                  tissue = tissue,
                  type = `Tumour type`,
                  mad2 = `Mad2 levels`,
                  stage = `Early/ Late`)
    index$tissue[is.na(index$tissue)] = "other"

    # load Gene Ontology sets from Ensembl
    go = gs$go(dset="mmusculus_gene_ensembl", genes="ensembl_gene_id") %>%
        filter(ensembl_gene_id %in% rownames(counts)) %>%
        transmute(gene = ensembl_gene_id,
                  set = paste(id, name)) %>%
        unstack()
    size = sapply(go, length)
    go = go[size >= 5 & size <= 200]

    deg = de_genes(counts, ~ tissue)
    gsa = de_sets(deg, go)
}
