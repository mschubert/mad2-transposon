setwd("C:/Users/bjorn/Dropbox/RNA-seq/012 - Sophia (SBTs)/SBT_FINAL/")

library(ggplot2)

# load annotated DE table from file # contains cpm read counts and gene annotations
load("SBT_-1_1_final_table.RData")

# head(final_table)
#                hgnc_symbol chromosome_name start_position end_position    position     logFC   logCPM       LR       PValue          FDR
#ENSG00000075275      CELSR1              22       46360834     46537170  46449002.0 -7.163614 6.226853 57.65728 3.120064e-14 4.391077e-10
#ENSG00000225972    MTND1P23               1         629062       629433    629247.5 -6.894271 6.354573 56.50251 5.612675e-14 4.391077e-10
#ENSG00000164434       FABP7               6      122779475    122784074 122781774.5 -7.294962 7.754199 53.33308 2.815261e-13 1.468346e-09
#ENSG00000115884        SDC1               2       20200797     20225433  20213115.0 -7.405456 4.911198 47.96070 4.348488e-12 1.701020e-08
#ENSG00000126562        WNK4              17       42780678     42796936  42788807.0 -7.399091 5.031155 46.39633 9.659647e-12 3.022890e-08
#ENSG00000108846       ABCC3              17       50634777     50692252  50663514.5 -7.586513 6.872393 45.97253 1.199226e-11 3.127381e-08

#                     SBT SBT15003  SBT15004  SBT15008
#ENSG00000075275 293.0468 1.760595 1.5447983 2.8011364
#ENSG00000225972 318.8662 1.942726 3.2749724 2.8011364
#ENSG00000164434 846.9906 7.406642 6.4881529 2.2647486
#ENSG00000115884 117.8589 1.214204 0.3707516 0.4767892
#ENSG00000126562 128.1013 1.274914 0.7415032 0.2383946
#ENSG00000108846 461.0508 3.521191 0.3707516 3.2779255

# only the chromosome_name, position, and counts for the samples are relevant (counts are CPM)


# filter final table to only have chromosomes 1 to 22 and X, and refactor for purpose of plotting
chromosomes <- c(1:22, 'X')
final_table <- final_table[final_table$chromosome_name %in% chromosomes, ]
final_table$chromosome_name <- factor(final_table$chromosome_name, levels = chromosomes)

# set sample IDs; SBT15004 = euploid reference control
sbts <- c("SBT15003", "SBT15004", "SBT15008")

# set file name
fileName <- "__RNA-seq_CNV_SBT15004_as_REF.pdf"

for(sbt in sbts) {
  
  # get log2 ratios of sample over reference; 
  final_table$ratio <- log((final_table[, sbt]+1)/(final_table$SBT15004+1), 2)
  
  # get median ratios per chromosome - will be plotted as red line in per chromosome plot
  medianRatios <- vector("numeric", length = length(chromosomes))
  names(medianRatios) <- chromosomes
  for(chromosome in chromosomes) {
    ratiosChr <- final_table$ratio[final_table$chromosome_name == chromosome]
    medianRatios[chromosome] <- median(ratiosChr)
  }
  hline.data <- data.frame(z = medianRatios, 
                           chromosome_name = chromosomes)
  
  pdf(file = paste(sbt, "_genome_wide", fileName), width = 15, height = 12)
  g <- ggplot(final_table, aes(x = position/1e+06, y = ratio)) +
    geom_point(alpha = 0.5) +
    geom_hline(aes(yintercept = z), hline.data, col = "red") +
    facet_wrap(~chromosome_name) +
    theme_bw() +
    labs(x = "Genomic position (Mb)", y = "log2 ratio of sample read count over euploid reference", 
         title = paste(sbt, " vs SBT15004 - by chromosome", sep = ""))
  print(g)
  dev.off()
  
  
  pdf(file = paste(sbt, "_per_chromosome", fileName), width = 7, height = 5)
  for(chromosome in chromosomes) {
    g <- ggplot(final_table[final_table$chromosome_name %in% chromosome, ], aes(x = position/1e+06, y = ratio)) +
      geom_point(alpha = 0.5) +
      geom_hline(aes(yintercept = z), hline.data[hline.data$chromosome_name == chromosome, ], col = "red") +
      theme_bw() +
      labs(x = "Genomic position (Mb)", y = "log2 ratio of sample read count over euploid reference", 
           title = paste(sbt, "vs SBT15004 - Chromosome", chromosome, sep = " "))
    print(g)
  }
  dev.off()
  
}  