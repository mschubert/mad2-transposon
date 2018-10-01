library(DESeq2)
library(ggplot2)
library(dplyr)
library(tidygraph)
io = import('io')
sys = import('sys')
set = import('./de_sets')
vp = import('../link_cis_expr/cor_viper')
idmap = import('process/idmap')
util = import('./util')

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'gene expression RData', 'de.RData'),
    opt('o', 'outfile', 'results RData', 'de_Ets+Erg.RData'),
    opt('p', 'plotfile', 'pdf', 'de_Ets+Erg.pdf'),
    arg('sets', 'gene set .RData', arity='*',
        list.files("../data/genesets", "\\.RData", full.names=TRUE)))

#om = function(x) idmap$orthologue(x, from="external_gene_name", to="mgi_symbol")
#net = io$load(args$network) %>%
#    mutate(Regulator = unname(om(Regulator)),
#           Target = unname(om(Target))) %>%
#    na.omit()
#tf_net = filter(net, Target %in% Regulator)

dset = io$load(args$diff_expr)
#cis = dset$cis$samples %>%
#    filter(external_gene_name == args$ins) %>%
#    select(sample, ins=external_gene_name)
eset = dset$eset
idx = colData(eset) %>%
    as.data.frame() %>%
#    left_join(cis) %>%
    mutate(#ins = ifelse(is.na(ins), 0, 1),
           aneup0.3 = pmin(aneuploidy, 0.3),
           keep = type %in% c("T-cell", "Other") & !sample %in% c("401t", "403t", "612t", "631s", "477t"),
           group = factor(DESeq2::counts(eset, normalized=TRUE)["ENSMUSG00000032035",] > 700))
levels(idx$group) = c("Erg", "Ets1")
idx$group = paste(idx$group, idx$tissue, sep=":")
idx$group[!idx$keep] = NA
idx$group = relevel(factor(idx$group), "Ets1:spleen")
eset@colData = DataFrame(idx)
eset = eset[,idx$keep]
idx = idx[idx$keep,]

#expr = assay(eset) #FIXME: normalized??
#rownames(expr) = idmap$gene(rownames(expr),
#    to="external_gene_name", dset="mmusculus_gene_ensembl")

res = do_wald(eset, ~ group)
res$aneup0.3 = do_lrt(eset, ~ group + aneup0.3, ~ group)
res$`T+Ets_aneup` = do_wald(eset[,idx$group == "Ets1:thymus"], ~ aneup0.3)
res$`B+Ets_aneup` = do_wald(eset[,idx$group == "Ets1:spleen"], ~ aneup0.3)
res$`B+Erg_aneup` = do_wald(eset[,idx$group == "Erg:spleen"], ~ aneup0.3)

sets = io$load(args$sets) %>%
    setNames(tools::file_path_sans_ext(basename(args$sets))) %>%
    lapply(function(x) set$gset$filter(x, min=5, valid=na.omit(res$gene_name)))

pdf(args$plotfile)
for (rname in names(res)) {
    message(rname)
    print(util$plot_volcano(res[[rname]]) + ggtitle(rname))
    for (sname in names(sets)) {
        title = paste(rname, sname)
        message(title)
        print(set$plot_gset(res[[rname]], sets[[sname]]) + ggtitle(title))
    }
}
dev.off()
#print(util$plot_pcs(idx, dset$pca, 1, 2, hl=cis$sample))
#
#dviper = vp$diff_viper(expr, net, eset$ins * eset$aneup0.2)
#dcor = vp$diff_cor(expr, tf_net, eset$ins * eset$aneup0.2)
#print(vp$plot_subnet(dviper, dcor) + ggtitle("MI network"))
#
#print(util$plot_volcano(res) + ggtitle("gene"))
#for (name in names(sets))
#    print(set$plot_gset(res, sets[[name]]) + ggtitle(name))

save(res, file=args$outfile)
