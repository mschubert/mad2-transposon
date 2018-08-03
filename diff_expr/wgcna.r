library(dplyr)
library(WGCNA)
io = import('io')
sys = import('sys')
idmap = import('process/idmap')
gnet = import('tools/genenet')
plt = import('plot')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression RData', '../data/rnaseq/assemble.RData'),
#    opt('c', 'copies', 'gene copy matrix', '../ploidy_compare/gene_copies.RData'),
    opt('c', 'cis', 'sites per sample', '../cis_analysis/poisson.RData'),
    opt('a', 'aneup', 'aneuploidy score', '../ploidy_compare/analysis_set.RData'),
    opt('o', 'outfile', 'results RData', 'wgcna.RData'),
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

pdf(args$plotfile, 15, 15)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(expr, powerVector = powers, verbose = 5)
plot(sft$fitIndices$Power, -sign(sft$fitIndices$slope)*sft$fitIndices$SFT.R.sq,
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",
     main = paste("Scale independence"))

adjacency = adjacency(expr, power=6)
TOM = TOMsimilarity(adjacency)
geneTree = hclust(as.dist(1-TOM), method = "average")
dynamicMods = cutreeDynamic(dendro = geneTree, distM = 1-TOM,
    deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = 30)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

# merge similar modules
MEList = moduleEigengenes(expr, colors = dynamicColors)
METree = hclust(as.dist(1-cor(MEList$eigengenes)), method = "average")
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=0.1, col = "red")
merged = mergeCloseModules(expr, dynamicColors, cutHeight=0.1, verbose = 3)
plotDendroAndColors(geneTree, cbind(dynamicColors, merged$colors),
    c("Dynamic Tree Cut", "Merged dynamic"),
    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05)

plot_trait_cor = function(traits, title, pval, abs=FALSE) {
    moduleTraitCor = cor(MEs, traits, use = "p")
    if (abs) {
        moduleTraitCor = abs(moduleTraitCor)
        col = plt$brew$seq(1:50)
    } else
        col = plt$brew$div(1:50)
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(idx))
    mtp = moduleTraitPvalue[,narray::map(moduleTraitPvalue, along=1, min) < pval]
    mtc = moduleTraitCor[,colnames(mtp)]

    ord_cols = colnames(mtc)[hclust(dist(t(mtc)))$order]
    ord_rows = rownames(mtc)[hclust(dist(mtc))$order]
    mtp = mtp[ord_rows, ord_cols]
    mtc = mtc[ord_rows, ord_cols]

    textMatrix = paste(signif(mtc, 2), "\n(",signif(mtp, 1), ")", sep = "")
    dim(textMatrix) = dim(mtc)
    labeledHeatmap(Matrix = mtc, xLabels = colnames(mtc), yLabels = names(MEs),
        ySymbols = names(MEs), colorLabels = FALSE, colors = col,
        textMatrix = textMatrix, cex.text=0.5, main = title)
}
MEs = orderMEs(merged$newMEs)
plot_trait_cor(MEs, "modules", pval=0.01)
plot_trait_cor(traits_ins, "insertions", pval=0.01)
plot_trait_cor(traits_expr, "expression", pval=0.01)
plot_trait_cor(traits_expr, "expression (absolute)", pval=0.01, abs=TRUE)
print(gnet$plot_pcor_net(gnet$pcor(t(traits_expr)), fdr=0.3))

both = cbind(MEs, traits_expr)
plot_trait_cor(both, "modules+expr", pval=0.01)
plot_trait_cor(both, "modules+expr (absolute cor)", pval=0.01, abs=TRUE)
#plot_trait_cor(both, "highlight ins expr", pval=0.01, abs=TRUE)
#print(gnet$plot_bootstrapped_pcor(t(both), fdr=0.3))
print(gnet$plot_pcor_net(gnet$pcor(t(both)), pval=0.05, node_text=4, edge_text=2))

dev.off()

save(merged, MEList, METree, MEs, file=args$outfile)
