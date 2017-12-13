library(dplyr)
library(cowplot)
library(ggrepel)
b = import('base')
io = import('io')

# load inferred ploidy for all samples
aneup = io$load('../aneuploidy/ploidy_eT.RData')$aneuploidy
names(aneup) = paste0(b$grep("([0-9]+)", names(aneup)),
                      toupper(b$grep("(s|t|S|T)", names(aneup))))

# load Mad2 expression
expr = io$load('../data/rnaseq/assemble.RData')$expr
colnames(expr) = paste0(b$grep("([0-9]+)", colnames(expr)),
                        toupper(b$grep("(s|t|S|T)", colnames(expr))))
mad2 = expr['ENSMUSG00000029910',]

# plot sample-level aneuploidy vs mad2 level
narray::intersect(aneup, mad2)

samples = data.frame(sample = names(aneup),
               aneup = aneup,
               mad2 = mad2) %>%
    mutate(mad2_class = ifelse(mad2 > 7, "high", "low"),
           aneup_class = ifelse(aneup > 7, "high", "low"))

p = ggplot(samples, aes(x=aneup, y=mad2)) +
    geom_point(aes(shape=mad2_class, color=aneup_class), size=5) +
    geom_text_repel(aes(label=sample)) +
    ggtitle("Inferred aneuploidy vs. Mad2 expression")

save(samples, file="aneuploidy_mad2.RData")

pdf("aneuploidy_mad2.pdf")
print(p)
dev.off()
