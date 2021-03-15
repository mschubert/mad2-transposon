library(dplyr)
sys = import('sys')
gset = import('genesets')
util = import('./util')

name2genes = function(cn) {
    cur_coll = colls[[cn]]
    cur_sets = cfg[[cn]]
    cur_coll[cur_sets]
}

args = sys$cmd$parse(
    opt('c', 'config', 'yaml', 'set-detail_Mad2PB.yaml'),
    opt('d', 'de_obj', 'rds', 'de_Mad2PB.rds'),
    opt('k', 'key', 'key in config', 'DNArepair'),
    opt('p', 'plotfile', 'pdf', 'de_Mad2PB_set-detail/MycV1.pdf')
)

cfg = yaml::read_yaml(args$config)[[args$key]]
colls = gset$get_mouse(names(cfg), drop=FALSE)
genes = lapply(names(colls), name2genes) %>%
    unlist() %>% unique()

res = readRDS(args$de_obj) %>%
    lapply(. %>% filter(gene_name %in% genes))

pdf(args$plotfile)
for (rname in names(res))
    print(util$plot_volcano(res[[rname]]) + ggtitle(rname))
dev.off()
