sys = import('sys')
gdsc = import('data/gdsc')

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', 'plot.yaml'),
    opt('i', 'id', 'key in yaml', 'brca'),
    opt('p', 'plotfile', 'pdf', 'brca.pdf'))

cfg = yaml::read_yaml(args$config)
genes = cfg$genes
tissues = cfg[[args$id]]$tissues
clines = cfg[[args$id]]$clines

tissues = gdsc$tissues(tissues)
expr = gdsc$basal_expression()
expr = expr[genes ,intersect(colnames(expr), names(tissues))]
colnames(expr) = gdsc$cosmic$id2name(colnames(expr))
aneup = gdsc$aneuploidy() %>% select(cline=cell_line, aneuploidy)

edf = reshape2::melt(expr) %>%
    transmute(gene=factor(Var1, levels=genes), cline=as.character(Var2), expr=value) %>%
    mutate(label = ifelse(cline %in% clines, cline, NA)) %>%
    left_join(aneup)

p = ggplot(edf, aes(x=gene, y=expr)) +
    ggbeeswarm::geom_quasirandom(aes(alpha=is.na(label), size=aneuploidy)) +
    ggrepel::geom_text_repel(aes(label=label)) +
    guides(alpha=FALSE) +
    scale_alpha_manual(values=c(1, 0.1))

pdf("gdsc_stat_myc.pdf", 12, 8)
print(p)
dev.off()
