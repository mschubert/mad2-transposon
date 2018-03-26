# Load required packages and dependencies
library(edgeR)
library(RColorBrewer)
library(ggplot2)
library(gplots)
library(geneplotter)
library(xlsx)
library(biomaRt)
library(EDASeq)


setwd("C:/Users/bjorn/Dropbox/RNA-seq/017 - Transposon screen/")

# Load count table
in_data1 <- read.table("PB170929_raw_readcount_dedup_UMIs.txt", header = TRUE, row.names = "Gene")
sampleIDs <- c("S157", "T180", "S404", "T407", "S409", "S411",
               "S413", "S416", "S417", "T418", "S422", "S424",
               "S428", "S429", "S431", "S434", "S435", "S437",
               "T449", "T477", "S482", "S412", "S451", "S476",
               "S489L", "T452", "S462", "T184", "T401", "T419",
               "S442", "T443")
colnames(in_data1) <- sampleIDs

data <- in_data1
remove(in_data1)

# Define groups, contrasts, and report (y/n)
projectTitle <- "Transposon_general"

# original grouping
sampleType <- c("spleen", "thymus", "spleen", "thymus", "spleen", "spleen",
                "spleen", "spleen", "spleen", "thymus", "spleen", "spleen",
                "spleen", "spleen", "spleen",  "spleen", "spleen", "spleen",
                "thymus", "thymus", "spleen", "spleen", "spleen", "spleen",
                "spleen", "thymus", "spleen", "thymus", "thymus", "thymus",
                "spleen", "thymus")

report <- 1

# tcr analysis
# Convert Ensembl IDs into gene symbols
library(biomaRt)
listMarts(host = "www.ensembl.org")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")
geneNames <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "start_position", "end_position"),
                   filters = "ensembl_gene_id",
                   values = rownames(data),
                   mart = mart)

# isolate the Trbv genes
trbv <- geneNames[grepl(pattern = "Trbv", x = geneNames$mgi_symbol), ]

# get transcript lengths
geneLengths <- getGeneLengthAndGCContent(trbv$ensembl_gene_id, org = "mmu")
geneNames <- cbind(trbv, geneLengths)

sampleIDs <- grepl(colnames(data), pattern = "T")

trbvCpm <- data[trbv$ensembl_gene_id, sampleIDs]
#trbvCount <- data[trbv$ensembl_gene_id, ]
row.names(trbvCpm) <- trbv$mgi_symbol
trbvCpm <- rpkm(trbvCpm, gene.length = geneNames$length, lib.size = colSums(data[sampleIDs]))
trbvCpmRel <- sweep(trbvCpm, 2, STATS = colSums(trbvCpm), FUN = "/")
cols <- rainbow(n = nrow(trbvCpmRel))
cols <- cols[sample(1:length(cols), size = length(cols), replace = FALSE)]
barplot(trbvCpmRel, col = cols)
barplot(trbvCpm, col = cols, cex.names = 0.5)

# isolate Bcr genes
# isolate the Trbv genes
brbv <- geneNames[grepl(pattern = "Brbv", x = geneNames$mgi_symbol), ]

# get transcript lengths
geneLengths <- getGeneLengthAndGCContent(trbv$ensembl_gene_id, org = "mmu")
geneNames <- cbind(trbv, geneLengths)

trbvCpm <- data[trbv$ensembl_gene_id, ]
#trbvCount <- data[trbv$ensembl_gene_id, ]
row.names(trbvCpm) <- trbv$mgi_symbol
trbvCpm <- rpkm(trbvCpm, gene.length = geneNames$length, lib.size = colSums(data))
trbvCpmRel <- sweep(trbvCpm, 2, STATS = colSums(trbvCpm), FUN = "/")
cols <- rainbow(n = nrow(trbvCpmRel))
cols <- cols[sample(1:length(cols), size = length(cols), replace = FALSE)]
barplot(trbvCpmRel, col = cols)
barplot(trbvCpm, col = cols)

# make copy of data as a back-up in case of accidents
#data.full <- data

# Making a metadata data frame
meta <- data.frame(
  row.names = colnames(data),
  condition = sampleType)
meta$condition <- relevel(meta$condition, ref = "spleen")

# Create a design matrix to specify the factors that are expected to affect expression levels
design <- model.matrix( ~ 0 + condition, meta)
contrasts <- list(c(-1, 1))



### \/ \/ \/ CHECK WHETHER CONTRASTS AND GROUPS ARE COHERENT, THEN AUTOMATE RUN FROM HERE \/ \/ \/ ###



# Create a DGEList object
d <- DGEList(counts = data, group = meta$condition)

# number of genes before filtering: 31294
dim(d)

d.full <- d # d.full acts as a back-up for d in case something goes wrong.
head(d$counts)
d$samples

finalGenes <- rowSums(cpm(d) > 1) > 2
d <- d.full[finalGenes, ]

# number of genes after filtering: 16291
dim(d)

# After filtering it is a good idea to reset the library size.

d$samples$lib.size <- colSums(d$counts)
d$samples # Should yield the same values as the total gene counts per sample


# Estimate normalization factors
d <- calcNormFactors(d, method = "TMM")
d$samples

# Generate a heatmap of samples, based on their similarity
norm_factors <- calcNormFactors(d, method="TMM")
pre_fpm <- sweep(data, 2, norm_factors$samples[, 2], '/')
fpm <- as.matrix(sweep(pre_fpm, 2, 1e6, '*'))

hmcols <- colorRampPalette(c("blue", "white", "red"))(n=255)
heatmap.2(fpm[1:1000, ], labRow = "", trace = "none", col = hmcols, scale = "row")


reducedFpm <- fpm[c("ENSMUSG00000029910",
                    "ENSMUSG00000059552",
                    "ENSMUSG00000024151",
                    "ENSMUSG00000013663"), ]
geneLengths <- getGeneLengthAndGCContent(id = row.names(reducedFpm), org = "mmu")

selectedGeneDataFrame <- data.frame(gene.id = row.names(reducedFpm),
                                    gene.symbol = c("Mad2l1", "Trp53", "Msh2", "Pten"),
                                    gene.length = geneLengths[, 1])
row.names(selectedGeneDataFrame) <- NULL
selectedGeneDataFrame <- cbind(selectedGeneDataFrame, reducedFpm)


reducedD <- d$counts[c("ENSMUSG00000029910",
                       "ENSMUSG00000059552",
                       "ENSMUSG00000024151",
                       "ENSMUSG00000013663"), ]
rpkmTable <- rpkm(x = cpm(reducedD, lib.size = d$samples$lib.size),
                  normalized.lib.sizes = TRUE,
                  lib.sizes = d$samples$lib.size,
                  gene.length = selectedGeneDataFrame$gene.length, log = FALSE)

par(mfrow = c(4, 1))
barplot(rpkmTable[1, ], col = "darkred", main = "Mad2l1")
barplot(rpkmTable[2, ], col = "darkred", main = "Trp53")
barplot(rpkmTable[3, ], col = "darkred", main = "Msh2")
barplot(rpkmTable[4, ], col = "darkred", main = "Pten")
par(mfrow = c(1, 1))


# PCA et al
pdf(file = paste0(projectTitle, "_report.pdf"), width = 12, height = 8)
heatmap(fpm[1:1000, ], labRow = "", Colv=FALSE, scale="row", rowlab = "", legend = T)

# Inspect relationship between samples using a multidimensional scaling plot of BCV and logFC
mdsCol <- brewer.pal(n = length(unique(sampleType)), name = "Set2")

par(mfrow = c(2, 3))
plotMDS(d, labels = colnames(data), col = mdsCol[factor = meta$condition],
        method = "bcv", main = "MDS: BCV, dim 1 v 2", dim.plot = c(1, 2))
plotMDS(d, labels = colnames(data), col = mdsCol[factor = meta$condition],
        method = "bcv", main = "MDS: BCV, dim 1 v 3", dim.plot = c(1, 3))
plotMDS(d, labels = colnames(data), col = mdsCol[factor = meta$condition],
        method = "bcv", main = "MDS: BCV, dim 2 v 3", dim.plot = c(2, 3))

plotMDS(d, labels = colnames(data), col = mdsCol[factor = meta$condition],
        method = "logFC", main = "MDS: logFC, dim 1 v 2", dim.plot = c(1, 2))
plotMDS(d, labels = colnames(data), col = mdsCol[factor = meta$condition],
        method = "logFC", main = "MDS: logFC, dim 1 v 3", dim.plot = c(1, 3))
plotMDS(d, labels = colnames(data), col = mdsCol[factor = meta$condition],
        method = "logFC", main = "MDS: logFC, dim 2 v 3", dim.plot = c(2, 3))
par(mfrow = c(1, 1))
dev.off()

# pca
rld <- rlog()


### COMPLEX DESIGN - using GLMs to correct for multiple factors ###


# Estimate the dispersion values, relative to the design matrix using Cox-Reid-adjusted likelihood
d <- estimateGLMCommonDisp(d, design)
d <- estimateGLMTrendedDisp(d, design)
d <- estimateGLMTagwiseDisp(d, design)


# Generate a visual representation of the mean-variance relationship
plotMeanVar(d, show.tagwise.vars = TRUE, NBline = TRUE, main = "Mean-variance relationship")
plotBCV(d, main = "Biological coefficient of variation vs mean log CPM")

# Given the design matrix and dispersion estimates, fit a GLM to each feature
fit_trended <- glmFit(d, design = design, dispersion = d$trended.dispersion)
fit_tagwise <- glmFit(d, design = design, dispersion = d$tagwise.dispersion)

fit <- fit_tagwise

# merging function, to be used in the loop
mergeAllData <- function(geneNames, ttTable, nc) {

  rn <- row.names(ttTable)
  merge1 <- merge(geneNames, ttTable, by = "row.names")
  row.names(merge1) <- merge1$Row.names
  merge1 <- merge1[, -1]

  merge1 <- merge(merge1, nc[rn, order(meta$condition)], by = "row.names")
  row.names(merge1) <- merge1$Row.names
  merge1 <- merge1[, -1]

  merge1 <- merge1[order(merge1$FDR, decreasing = FALSE), ]
  return(merge1)

}

options(java.parameters = "-Xmx8000m")

for(i in 1:length(contrasts)) {

  contrast <- contrasts[[i]]

  # Perform a likelihood ratio test, specifying the difference of interest (i.e. the colum number 'design')
  de <- glmLRT(fit, contrast = contrast)

  # Present a tabular summary of the differential expression statistics
  tt <- topTags(de, n = nrow(d))

  head(tt$table)

  # Inspect the depth-adjusted reads per million for several top differentially expressed genes
  nc <- cpm(d, normalized.lib.sizes = TRUE)
  rn <- rownames(tt$table)
  head(nc[rn, order(meta$condition)], 6)

  # Convert Ensembl IDs into gene symbols
  library(biomaRt)
  listMarts(host = "www.ensembl.org")
  mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")
  geneNames <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "start_position", "end_position"),
                     filters = "ensembl_gene_id",
                     values = rownames(tt$table),
                     mart = mart)

  position <- (geneNames$start + geneNames$end)/2
  geneNames <- cbind(geneNames, position)
  geneNames <- geneNames[!duplicated(geneNames$ensembl_gene_id), ]
  row.names(geneNames) <- geneNames$ensembl_gene_id
  geneNames <- geneNames[, -1]

  fileName <- paste(projectTitle, as.character(paste(contrast, collapse = "_")), "logFC_1_final_table", sep = "_")

  final_table <- mergeAllData(geneNames, tt$table, nc)
  save(final_table, file = paste(fileName, ".RData", sep = ""))
  write.table(final_table, file = paste(fileName, ".txt", sep = ""), quote = FALSE, dec = ",", sep = "\t", row.names = TRUE, col.names = TRUE)
  #write.xlsx2(final_table, file = paste(fileName, ".xlsx", sep = ""), dec = ",", sep = "\t", row.names = TRUE, col.names = TRUE)
  write.xlsx2(final_table[final_table$FDR < 0.05, ], file = paste(fileName, "sig_only.xlsx", sep = ""), dec = ",", sep = "\t", row.names = TRUE, col.names = TRUE)

}


# Plot a group correlation matrix
correlation <- cor(predFC(d, design, prior.count = 1, dispersion = 0.05))
print(correlation)
hmcol = colorRampPalette(c("darkblue", "white"))(n=255)
heatmap.2(correlation[-9, -9], trace = "none", col = hmcol, margins = c(10, 10))

dev.off()


# get gene length based on cpm
library(EDASeq)
geneLengths <- getGeneLengthAndGCContent(id = row.names(d$counts)[1:100], org = "mmu")


tableFiles <- list.files("raw comparison files/", pattern = ".RData", full = TRUE)
load(tableFiles[1])
dataMat <- final_table[final_table$FDR < 0.05, -c(2:10)]
rownames(dataMat) <- dataMat[, 1]
dataMat <- as.matrix(dataMat[, -1])
heatmap.2(dataMat, col = hmcols, trace = "none", scale = "row", labRow = "")



# k-means cluster


library(amap)
load("T-ALL_general_v2_-1_0_1_0_0_0_0_logFC_1_final_table.RData")
data.xpr <- final_table[, 11:33]
numClust <- 20
clusters = Kmeans(x = data.xpr, centers=numClust, method="pearson", iter.max = 100)

clustersNames <- as.numeric(clusters$cluster)
scaledMat <- t(apply(data.xpr, 1, function(x) {as.vector(scale(x))}))
colnames(scaledMat) <- colnames(data.xpr)

geneClustList <- list()
sampleOrder <- colnames(data.xpr)
for(i in 1:numClust) {
  geneClustList[[paste0("cluster_", i)]] <- scaledMat[clustersNames == i, sampleOrder]
}

hmcol <- colorRampPalette(c("blue", "white", "red"))(n=255)
pdf(file = paste0(numClust, "_k-means_clusters_170919.pdf"))
for(i in names(geneClustList)) {
  heatmap.2(geneClustList[[i]], scale = "none", col = hmcol, Colv = NA, dendogram = "row",
            trace = "none", density.info = "none", labRow = "", main = i)
}
dev.off()

groups_supervised <- list(c("T372", "T372"),
                          c("T317", "T318", "T503"),
                          c("eT_p0", "eT_p2"),
                          c("T559", "T560", "T571"),
                          c("T363", "T496", "T548", "T553", "T554"),
                          c("T165", "T179", "T187", "T223", "T230", "T231", "T330", "T362"))

groupMeans <- matrix(NA, nrow = nrow(data.xpr), ncol = length(groups_supervised))
for(i in 1:length(groups_supervised)) {
  for(j in 1:nrow(scaledMat)) {
    groupMeans[j, i] <- mean(scaledMat[j, groups_supervised[[i]]])
  }
}
colnames(groupMeans) <- c("wt", "aneuploid", "eT", "lo", "mid", "T-ALL")
rownames(groupMeans) <- rownames(data.xpr)

groupMeanClustered <- list()
for(i in 1:numClust) {
  groupMeanClustered[[paste0("cluster_", i)]] <- groupMeans[clustersNames == i, ]
}

pdf(file = paste0(numClust, "_k-means_cluster_profiles_new_v2.pdf"), height = 16, width = 15)
par(mfrow = c(4, 5))
for(i in 1:numClust) {
  plotMat <- groupMeanClustered[[i]]
  plot(plotMat[1, ], type = "l", main = paste0("cluster ", i), ylim = c(-3.5, 3.5), xaxt = "n",
       ylab = "normalized expression", xlab = "group", lwd = 0.1, col=alpha(rgb(0,0,0), 0.3))
  for(j in 2:nrow(plotMat)) {
    lines(plotMat[j, ], lwd = 0.1, col=alpha(rgb(0,0,0), 0.3), xaxt = "n")
  }
  lines(colMeans(plotMat), lwd = 2, col = "red", xaxt = "n")
  axis(side = 1, labels = colnames(groupMeans), at = 1:length(groups_supervised))
}
dev.off()

for(name in names(groupMeanClustered)) {
  write.xlsx(x = groupMeanClustered[[name]], append = TRUE, row.names = TRUE, col.names = TRUE,
             file = paste0(numClust, "cluster_gene_list_new.xlsx"), sheetName = name)
  write.table(x = groupMeanClustered[[name]], row.names = TRUE, col.names = TRUE, quote = FALSE,
              sep = "\t", file = paste0(numClust, "_", name, "_cluster_gene_list.txt"))
}

heatmap.2(groupMeanClustered[[14]][test$userid, ], scale = "none", col = hmcol, Colv = NA, dendogram = "row",
          trace = "none", density.info = "none", labRow = test$Gene.Symbol)
