library(ggplot2)
library(dplyr)
b = import('base')
io = import('io')
idmap = import('process/idmap')

# load the gene expression data + sample types
#aneuploidy = io$load('../aneuploidy/ploidy_from_rnaseq.RData')$aneuploidy
dset = io$load('../data/rnaseq/assemble.RData')
expr = dset$expr[rowSums(dset$counts) >= 10,]
index = dset$idx %>%
    transmute(id = id,
              tissue = tissue,
              type = `Tumour type`,
              mad2 = `Mad2 levels`,
              stage = `Early/ Late`) #%>%
#    full_join(data.frame(id=names(aneuploidy), aneuploidy=unname(aneuploidy)), by = "id")
index$tissue[is.na(index$tissue)] = "other" # should this be spleen?

edf = reshape2::melt(expr) %>%
    transmute(ensembl_gene_id = as.character(Var1),
              sample = as.character(Var2),
              expr = value,
              gene = idmap$gene(ensembl_gene_id, to="external_gene_name",
                                dset="mmusculus_gene_ensembl")) %>%
    mutate(label = ifelse(gene %in% c("Erg", "Ets1"), gene, NA))

pdf("ets1_erg.pdf", 12, 5)
ggplot(edf, aes(x=sample, y=expr)) +
    geom_violin() +
    geom_text(aes(label=label)) +
    theme_bw() +
	labs(title="Expression level of Ets1/Erg compared to other genes",
		 subtitle="Not normalized by gene length")
dev.off()
