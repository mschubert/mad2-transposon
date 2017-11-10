# Load required packages and dependencies
library(edgeR)
library(biomaRt)
library(EDASeq)

# there should be an object called 'data' containing the raw (unfiltered) read counts


# load all transcripts and raw counts
load("all_transcripts_sampleOrdered.RData")

# head(data)
#                   T372 T373 T317 T318 T503 T559 T571 T363 T560 T496 T554 T362 T553 T165 T179 T187 T223 T230 T231 T330 T548 eT_p0 eT_p2
#ENSMUSG00000000001 1533 2479 1688 2585 3380 2703 2715 3879 1700 4161 1471 3068 4896 3103 4377 2501 2949 2927 3356 3020 1240  6191  4740
#ENSMUSG00000000003    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0    0     0     0
#ENSMUSG00000000028  230  421  453  643  782 1193 1461 1493  615 1778  913 1487 2654 1441  534 1142  894 1283 1496 1488  924  1334  1922
#ENSMUSG00000000031   33   50  109  161   81   73 2471   64 1395 1484  185  298   37  304  134   27 1440  135  531 1190    0    18  1642
#ENSMUSG00000000037   19   47   32   28   50   18    4    2    2    0   12    1    5    0   37    1    5    5    1    3    0     1     0
#ENSMUSG00000000049    0    0    0    0    1    0    0    0    0    0    0    0    1    0    0    1    0    1    0    0    0     0     0

# annotate the genes
listMarts(host = "www.ensembl.org")
mart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "www.ensembl.org")
geneNames <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "chromosome_name", "start_position", "end_position"),
                   filters = "ensembl_gene_id",
                   values = rownames(data),
                   mart = mart)

# isolate the Trbv genes
trbv <- geneNames[grepl(pattern = "Trbv", x = geneNames$mgi_symbol), ]

# alternative patterns for different Tcr loci: Trbj, Trg, Trd

# get transcript lengths
geneLengths <- getGeneLengthAndGCContent(trbv$ensembl_gene_id, org = "mmu") # select mmusculus as genome upon prompt
geneNames <- cbind(trbv, geneLengths)

# convert trbv genes to rpkm
trbvCpm <- data[trbv$ensembl_gene_id, ]
row.names(trbvCpm) <- trbv$mgi_symbol
trbvCpm <- rpkm(trbvCpm, gene.length = geneNames$length, lib.size = colSums(data))

# get relative contributions and plot
trbvCpmRel <- sweep(trbvCpm, 2, STATS = colSums(trbvCpm), FUN = "/")
cols <- rainbow(n = nrow(trbvCpmRel))
cols <- cols[sample(1:length(cols), size = length(cols), replace = FALSE)]
barplot(trbvCpmRel, col = cols)

