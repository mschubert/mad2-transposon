library(ggplot2)
library(ggrepel)
io = import('io')

et = io$load('ploidy_eT.RData')$dfs
tall = io$load('ploidy_cancer_only.RData')$ctl3 %>%
    transmute(chr = seqnames,
              sample = sample,
              rna2 = ploidy,
              wgs2 = ref_ploidy)

conservative_deviation = function(x, y) {
    x = c(x,y)
    x = x - 2
    minx = min(x)
    maxx = max(x)
    if (minx * maxx <= 0)
        2
    else if (minx < 0)
        maxx + 2
    else if (minx >= 0)
        minx + 2
    else
        stop("error")
}

both = dplyr::inner_join(et, tall, by=c('sample', 'chr'))
both$rna_both = purrr::map2_dbl(both$rna, both$rna2, conservative_deviation)
stopifnot(both$wgs == both$wgs2)
both$wgs2 = NULL

mod = lm(rna_both ~ wgs, data=both)
p = ggplot(both, aes(x=wgs, y=rna_both)) +
    geom_smooth(method="lm") +
    geom_point() +
    geom_text_repel(aes(label=label), size=3) +
    labs(x = "Ploidy from single-cell WGS data + AneuFinder",
         y = "Ploidy inferred from RNA-seq",
         title = sprintf("RNA-seq for ploidy inference (p=%.2g, r^2=%.2f)",
             mod %>% broom::tidy() %>% filter(term == "wgs") %>% pull(p.value),
             mod %>% broom::glance() %>% pull(r.squared)),
         subtitle = "assuming total DNA=const, diff euploid lib as reference")

pdf("both_conservative.pdf")
print(p)
dev.off()
