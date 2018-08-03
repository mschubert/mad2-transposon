library(dplyr)
library(WGCNA)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')
gnet = import('tools/genenet')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression RData', '../data/rnaseq/assemble.RData'),
#    opt('c', 'copies', 'gene copy matrix', '../ploidy_compare/gene_copies.RData'),
    opt('c', 'cis', 'sites per sample', '../cis_analysis/poisson.RData'),
    opt('a', 'aneup', 'aneuploidy score', '../ploidy_compare/analysis_set.RData'),
    opt('o', 'outfile', 'results RData', 'aneup_de.RData'),
    opt('p', 'plotfile', 'pdf', 'wgcna.pdf'))

aneup = io$load(args$aneup)
ins = io$load(args$cis)
is_cis = ins$result %>% filter(adj.p < 1e-3) %>% pull(external_gene_name) %>% unique()
is_cis_strict = ins$result %>% filter(adj.p < 1e-10) %>% pull(external_gene_name) %>% unique()
is_cis_strict = is_cis_strict[!grepl("^Gm", is_cis_strict)]
cis = ins$samples %>%
    filter(external_gene_name %in% is_cis) %>%
    narray::construct(n_ins ~ sample + external_gene_name, fill=0)
cis = (!is.na(cis) & cis!=0) + 0

eset = io$load(args$expr) # filter on read count, variance?
rownames(eset$expr) = eset$gene
expr = t(eset$expr[rowMeans(eset$counts >= 10) >= 0.1 | eset$gene %in% is_cis,])
idx = eset$idx %>% select(sample, tissue, type)
narray::intersect(expr, idx$sample, aneup$sample, cis, along=1) # loses 7 samples
pca = prcomp(t(expr), scale=TRUE)
traits_ins = cbind(pca$rotation[,1:4], aneup=aneup$aneup, cis)
traits_expr = cbind(aneup=aneup$aneup, expr[,is_cis_strict])

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(expr, powerVector = powers, verbose = 5)

pdf(args$plotfile, 10, 7)

plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope)*sft$fitIndices$SFT.R.sq,
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     main = paste("Scale independence"))

net = blockwiseModules(expr, power = sft$powerEstimate,
    TOMType = "unsigned", minModuleSize = 30,
    reassignThreshold = 0, mergeCutHeight = 0.25,
    numericLabels = TRUE, pamRespectsDendro = FALSE,
    verbose = 3)

mergedColors = labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
    "Module colors",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05)

MEs0 = moduleEigengenes(expr, net$colors)$eigengenes
MEs = orderMEs(MEs0)

plot_trait_cor = function(traits, title, pval) {
    moduleTraitCor = cor(MEs, traits, use = "p")
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(idx))

    mtp = moduleTraitPvalue[,narray::map(moduleTraitPvalue, along=1, min) < pval]
    mtc = moduleTraitCor[,colnames(mtp)]

    ord_cols = colnames(mtc)[hclust(dist(t(mtc)))$order]
    ord_rows = rownames(mtc)[hclust(dist(mtc))$order]
    mtp = mtp[ord_rows, ord_cols]
    mtc = mtc[ord_rows, ord_cols]

    textMatrix = paste(signif(mtc, 2), "\n(",signif(mtp, 1), ")", sep = "")
    dim(textMatrix) = dim(mtc)
    labeledHeatmap(Matrix = mtc,
        xLabels = colnames(mtc),
        yLabels = names(MEs),
        ySymbols = names(MEs),
        colorLabels = FALSE,
        colors = greenWhiteRed(50),
        textMatrix = textMatrix,
        setStdMargins = FALSE,
        cex.text = 0.5,
        zlim = c(-1,1),
        main = title)
}
plot_trait_cor(MEs, "modules", pval=0.01)
plot_trait_cor(traits_ins, "insertions", pval=0.01)
plot_trait_cor(traits_expr, "expression", pval=0.01)
print(gnet$plot_pcor_net(gnet$pcor(t(traits_expr)), fdr=0.2))

both = cbind(MEs, traits_expr)
plot_trait_cor(both, "modules+expr", pval=0.01)
print(gnet$plot_pcor_net(gnet$pcor(t(both)), fdr=0.1))

dev.off()
