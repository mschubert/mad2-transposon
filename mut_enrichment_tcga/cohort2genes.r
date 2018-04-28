library(dplyr)
b = import('base')
sys = import('sys')
plt = import('plot')
tcga = import('data/tcga')

do_fit = function(gene) {
    mdata = mut %>%
        filter(Hugo_Symbol == gene) %>%
        select(Sample) %>%
        mutate(mut = 1) %>%
        right_join(sdata, by="Sample") %>%
        mutate(mut = ifelse(is.na(mut), 0, mut))

    mod = glm(aneuploidy ~ n_mut + mut, family=Gamma(), data=mdata) %>%
        broom::tidy() %>%
        filter(term == "mut") %>%
        select(-term)
}

args = sys$cmd$parse(
    opt('c', 'cohort', 'TCGA cohort', 'LUAD'),
    opt('o', 'outfile', 'file to save model to', 'LUAD.RData'),
    opt('p', 'plotfile', 'pdf for plotting', 'LUAD.pdf'))

aneup = tcga$aneuploidy(cohort=args$cohort)
mut = tcga$mutations() %>%
    filter(Study == args$cohort) %>%
    rename(Sample = Tumor_Sample_Barcode) %>%
    tcga$filter(along="Sample", primary=TRUE, cancer=TRUE)

sdata = mut %>%
    group_by(Sample) %>%
    summarize(n_mut = n()) %>%
    inner_join(aneup)

result = mut %>%
    group_by(Hugo_Symbol) %>%
    summarize(n = n()) %>%
    filter(n >= 10) %>%
    mutate(result = purrr::map(Hugo_Symbol, do_fit)) %>%
    tidyr::unnest() %>%
    mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
    arrange(p.value)

p = result %>%
    mutate(size = n,
           label = Hugo_Symbol) %>%
    plt$p_effect("adj.p") %>%
    plt$volcano(base.size=0.5, repel=TRUE)

x = mut %>%
    select(Sample, Hugo_Symbol) %>%
    mutate(mut = TRUE) %>%
    tidyr::complete(Sample, Hugo_Symbol, fill=list(mut=FALSE)) %>%
    inner_join(sdata %>% select(Sample, aneuploidy)) %>%
    filter(Hugo_Symbol %in% head(result$Hugo_Symbol, 10))
p2 = ggplot(x, aes(x=Hugo_Symbol, y=aneuploidy)) +
    ggbeeswarm::geom_quasirandom(aes(fill=mut), shape=21, alpha=0.5) +
    theme(axis.text.x = element_text(angle=45, hjust=1))

pdf(args$plotfile)
print(p)
print(p2)
dev.off()

save(result, file=args$outfile)
