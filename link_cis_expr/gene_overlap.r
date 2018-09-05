library(patchwork)
io = import('io')
sys = import('sys')
plt = import('plot')

args = sys$cmd$parse(
    opt('s', 'select', 'yaml which sets to use', './interesting_sets.yaml'),
    opt('p', 'plotfile', 'pdf', 'gene_overlap.pdf'),
    arg('genesets', 'RData files', arity='*',
        list.files("../data/genesets", "\\.RData$", full.names=TRUE))
)

select = io$read_yaml(args$select)$expr_sets
sets = io$load(args$genesets)
sets = lapply(names(select), function(s) sets[[s]][select[[s]]]) %>%
    unlist(recursive=FALSE, use.names=TRUE)

jaccard = expand.grid(set1 = names(sets), set2=names(sets),
        stringsAsFactors=FALSE) %>%
    mutate(j = purrr::map2_dbl(set1, set2,
        function(a, b) length(intersect(sets[[a]],sets[[b]])) /
            length(union(sets[[a]],sets[[b]])))) %>%
    plt$cluster(j ~ set1 + set2)

p1 = data.frame(set = factor(names(sets), levels=levels(jaccard$set1))) %>%
    mutate(n = purrr::map_dbl(set, function(s) length(sets[[s]]))) %>%
    ggplot(aes(x=set, y=n)) +
        geom_col() +
        plt$theme$no_gx()

p2 = plt$matrix(jaccard, j ~ set1 + set2, palette="Blues") +
    coord_fixed()

pdf(args$plotfile, 10.5, 10)
print(p1 + p2 + plot_layout(ncol=1, heights=c(1,5)))
dev.off()
