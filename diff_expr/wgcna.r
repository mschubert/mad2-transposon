library(dplyr)
library(WGCNA)
io = import('io')
sys = import('sys')

args = sys$cmd$parse(
    opt('e', 'expr', 'gene expression RData', '../data/rnaseq/assemble.RData'),
#    opt('c', 'copies', 'gene copy matrix', '../ploidy_compare/gene_copies.RData'),
    opt('c', 'cis', 'sites per sample', '../cis_analysis/poisson.RData'),
    opt('a', 'aneup', 'aneuploidy score', '../ploidy_compare/analysis_set.RData'),
    opt('o', 'outfile', 'results RData', 'aneup_de.RData'),
    opt('p', 'plotfile', 'pdf', 'aneup_de.pdf'))

aneup = io$load(args$aneup)
cis = io$load(args$cis)
is_cis = cis$result %>% filter(adj.p < 1e-3) %>% pull(external_gene_name)
traits = cis$samples %>%
    filter(external_gene_name %in% is_cis) %>%
    narray::construct(n_ins ~ sample + external_gene_name, fill=0)
traits = (!is.na(traits) & traits!=0) + 0

eset = io$load(args$expr) # filter on read count, variance?
expr = t(eset$expr[rowMeans(eset$counts >= 10) >= 0.1,])
idx = eset$idx %>% transmute(sample=paste0(hist_nr, tissue), tissue, type)
narray::intersect(expr, idx$sample, traits, along=1) # loses 7 samples
pca = prcomp(t(expr), scale=TRUE)
traits = cbind(pca$rotation[,1:4],
    aneup=aneup$aneup[match(rownames(traits), aneup$sample)], traits)

powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(expr, powerVector = powers, verbose = 5)

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
moduleTraitCor = cor(MEs, traits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nrow(idx))

mtp = moduleTraitPvalue[,narray::map(moduleTraitPvalue, along=1, min) < 1e-2]
mtc = moduleTraitCor[,colnames(mtp)]

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
    main = paste("Module-trait relationships"))

dev.off()
