library(dplyr)
library(DESeq2)
library(survival)
library(survminer)
tcga = import('data/tcga')
idmap = import('process/idmap')

do_wald = function(eset, fml, ex=NULL, drop=TRUE) {
    design(eset) = fml
    res = DESeq2::estimateDispersions(eset) %>%
        DESeq2::nbinomWaldTest(maxit=1000)
    if (length(ex) == 0)
        ex = setdiff(DESeq2::resultsNames(res), "Intercept")
    else
        ex = grep(ex, DESeq2::resultsNames(res), value=TRUE)
    res = sapply(ex, extract_coef, res=res, simplify=FALSE)
    if (length(res) == 1 && drop)
        res = res[[1]]
    res
}
extract_coef = function(res, coef, type="apeglm") {
    DESeq2::lfcShrink(res, coef=coef, type=type) %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene_name") %>%
        as_tibble() %>%
        mutate(stat = log2FoldChange / lfcSE) %>%
        arrange(padj, pvalue)
}
over_dmso = function(genotype, time, treatment) {
    cur = eset[,eset$time %in% time & eset$genotype == genotype &
                eset$treatment %in% c(treatment, "dmso")]
    colData(cur) = droplevels(colData(cur))
    res = do_wald(cur, ~ treatment, ex="treatment")
}
over_wt = function(genotype, time, treatment) {
    cur = eset[,eset$time %in% time & eset$genotype %in% genotype &
                eset$treatment %in% treatment]
    colData(cur) = droplevels(colData(cur))
    res = do_wald(cur, ~ genotype, ex="genotype")
}

eset = DESeq2::estimateSizeFactors(readRDS("dset.rds")$eset)
res = list()
res$wt_rev48_over_dmso = over_dmso("wt", "48", "rev")
res$rev48_stat1_over_wt = over_wt(c("wt", "stat1"), "48", "rev")
old_de = readRDS("../../expr_stat1/diff_expr.rds")[
    c("wt_rev48_over_dmso", "rev48_stat1_over_wt")]
x = inner_join(res$wt_rev48_over_dmso %>% select(gene_name, log2FoldChange),
               old_de$wt_rev48_over_dmso %>% select(gene_name, log2FoldChange),
               by="gene_name")
y = inner_join(res$rev48_stat1_over_wt %>% select(gene_name, log2FoldChange),
               old_de$rev48_stat1_over_wt %>% select(gene_name, log2FoldChange),
               by="gene_name")
plot(x$log2FoldChange.x, x$log2FoldChange.y) # same
plot(y$log2FoldChange.x, y$log2FoldChange.y) # same(!!)

sets = lapply(res, . %>% filter(log2FoldChange > 0) %>% pull(gene_name) %>% head(70))
old_sets = readRDS("../genesets/human/stat1_ko.rds")[
    c("wt_rev48_over_dmso", "rev48_stat1_over_wt")]
length(intersect(sets$wt_rev48_over_dmso, old_sets$wt_rev48_over_dmso)) # 70 of 70
length(intersect(sets$rev48_stat1_over_wt, old_sets$rev48_stat1_over_wt)) # 70 of 70
expr = tcga$rna_seq("BRCA", trans="vst")
rownames(expr) = idmap$gene(rownames(expr), to="hgnc_symbol")
scores = t(GSVA::gsva(expr, sets)) #TODO: scores are not used yet in brca/dset after
old_gsva = t(readRDS("../gsva/tcga-brca/stat1_ko.rds"))
narray::intersect(dset, scores, old_gsva, along=1)
plot(scores[,"wt_rev48_over_dmso"], dset$wt_rev48_over_dmso) # ~same
plot(scores[,"rev48_stat1_over_wt"], dset$rev48_stat1_over_wt) # ~same
plot(scores[,"wt_rev48_over_dmso"], old_gsva[,"wt_rev48_over_dmso"]) # ~same
plot(scores[,"rev48_stat1_over_wt"], old_gsva[,"rev48_stat1_over_wt"]) # ~same
plot(dset$rev48_stat1_over_wt, old_gsva[,"rev48_stat1_over_wt"]) # same

brca = readRDS("../../tcga_myc/dset.rds")
brca$meta$vital_status = as.integer(brca$meta$vital_status) - 1
dset = cbind(brca$meta, as.data.frame(brca$dmat)) %>%
    mutate(OS_years = OS_time / 365,
           vital_status = ifelse(OS_years > 10, 0, vital_status),
           OS_years = pmin(OS_years, 10)) %>%
    mutate(iclass = case_when(
        wt_rev48_over_dmso<0 & rev48_stat1_over_wt>0 ~ "CIN_stat1ko",
        wt_rev48_over_dmso>0 ~ "CIN",
        wt_rev48_over_dmso<0 & rev48_stat1_over_wt<0 ~ "noCIN"
    ), iclass=relevel(factor(iclass), "noCIN"))

p53wt = dset %>% filter(p53_mut == 0)
m1 = coxph(Surv(OS_years, vital_status) ~ age_at_diagnosis + purity + iclass, data=p53wt)
m1p = broom::tidy(m1) %>% filter(term == "iclassCIN_stat1ko") %>% pull(p.value)
fit1 = survfit(Surv(OS_years, vital_status) ~ iclass, data=p53wt)
ggsurvplot(fit1, data=p53wt, xlim=c(0,10), break.time.by=2.5)$plot +
    ylim(c(0.25,1)) +
    labs(subtitle = sprintf("p53 wt (n=%i) p=%.2g", sum(fit1$n), m1p))
