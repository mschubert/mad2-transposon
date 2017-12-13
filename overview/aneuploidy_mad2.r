library(dplyr)
library(cowplot)
library(ggrepel)
b = import('base')
io = import('io')

# load inferred ploidy for all samples
aneup = io$load('../aneuploidy/ploidy_eT.RData')$aneuploidy
names(aneup) = paste0(toupper(b$grep("(s|t|S|T)", names(aneup))),
                      b$grep("([0-9]+)", names(aneup)))

# load Mad2 expression
expr = io$load('../data/rnaseq/assemble.RData')$expr
colnames(expr) = paste0(toupper(b$grep("(s|t|S|T)", colnames(expr))),
                        b$grep("([0-9]+)", colnames(expr)))
mad2 = expr['ENSMUSG00000029910',]

# plot sample-level aneuploidy vs mad2 level
narray::intersect(aneup, mad2)

p = data.frame(sample = names(aneup),
               aneup = aneup,
               mad2 = mad2) %>%
    ggplot(aes(x=aneup, y=mad2)) +
        geom_point() +
        geom_text_repel(aes(label=sample)) +
        ggtitle("Inferred aneuploidy vs. Mad2 expression")

pdf("aneuploidy_mad2.pdf")
print(p)
dev.off()
