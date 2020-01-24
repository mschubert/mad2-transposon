sys = import('sys')
gdsc = import('data/gdsc')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', 'plot.yaml'),
    opt('i', 'id', 'key in yaml', 'brca'),
    opt('p', 'plotfile', 'pdf', 'brca.pdf'))

load_set = function(yset, hl) {
    x = readRDS(file.path("../data/geneset_scores/gdsc/", paste0(yset, ".rds")))
    re = x[hl[[yset]], ,drop=FALSE]
}

cfg = yaml::read_yaml(args$config)
hl = cfg$highlight
genes = hl$genes
tissues = cfg[[args$id]]$tissues
clines = cfg[[args$id]]$clines
sets = lapply(setdiff(names(hl), "genes"), load_set, hl=hl) %>%
    narray::stack(along=2)

tissues = gdsc$tissues(tissues)
expr = gdsc$basal_expression()
expr = expr[intersect(rownames(expr), genes),
            intersect(colnames(expr), names(tissues))]
eset = sets[,intersect(colnames(expr), names(tissues))]
#eset = narray::map(eset, along=2, scale)
expr = narray::stack(list(expr, eset), along=1)
colnames(expr) = gdsc$cosmic$id2name(colnames(expr))
aneup = gdsc$aneuploidy() %>% select(cline=cell_line, aneuploidy)

edf = reshape2::melt(expr) %>%
    transmute(gene=Var1, cline=as.character(Var2), expr=value,
              is_gene=ifelse(gene %in% genes, "gene", "geneset")) %>%
    mutate(label = ifelse(cline %in% clines, cline, NA)) %>%
    left_join(aneup)

p = ggplot(edf, aes(x=gene, y=expr)) +
    ggbeeswarm::geom_quasirandom(aes(alpha=is.na(label), size=aneuploidy)) +
    ggrepel::geom_text_repel(aes(label=label)) +
    guides(alpha=FALSE) +
    facet_wrap(~ is_gene, ncol=1, scales="free") +
    scale_alpha_manual(values=c(1, 0.1)) +
    theme(axis.text.x = element_text(angle=10, hjust=1))

pdf(args$plotfile, 14, 12)
print(p)
dev.off()
