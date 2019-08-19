library(dplyr)
sys = import('sys')
import('../../cgas/cor_structure/genenet', attach=TRUE)

args = sys$cmd$parse(
    opt('i', 'infile', 'rds', 'scores.rds'),
    opt('p', 'plotfile', 'pdf', 'genenet.pdf'))

data = readRDS(args$infile)

mat = t(data.matrix(data[,! colnames(data) %in% c("sample", "cohort")]))
colnames(mat) = data$sample
groups = narray::mask(data$cohort) + 0
tmat = narray::split(mat, along=2, subsets=data$cohort)
mat2 = rbind(mat, t(groups))

pdf(args$plotfile, 20, 15)
plot_cor_matrix(t(mat2), text_color=NULL)
plot_cor_matrix(t(mat), text_color=NULL)

pcor(t(mat)) %>% plot_pcor_net(node_size=4, edge_size=2.5)
plot_bootstrapped_pcor(t(mat), node_size=4)
pcor(t(mat2)) %>% plot_pcor_net(node_size=4, edge_size=2.5)
plot_bootstrapped_pcor(t(mat2), node_size=4)

for (i in seq_along(tmat)) {
    name = names(tmat)[i]
    plot_cor_matrix(t(tmat[[i]]), title=name, text_color=NULL)
    plot_bootstrapped_pcor(t(tmat[[i]]), fdr=0.3, node_size=4, title=name)
}
dev.off()
