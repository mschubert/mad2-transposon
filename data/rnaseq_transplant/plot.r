library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')

plot_one = function(title, comp) {
    plots = mapply(function(c, n) plt$volcano(c, label_top=30, repel=TRUE) + ggtitle(n),
        c=comp, n=names(comp), SIMPLIFY=FALSE)
    plt$text(title) / wrap_plots(plots, nrow=1) + plot_layout(heights=c(1,20))
}

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'rds', 'diff_expr.rds'),
    opt('p', 'plotfile', 'pdf', 'diff_expr.pdf')
)

res = readRDS(args$diff_expr)

pdf(18, 8, file=args$plotfile)
for (i in seq_along(res))
    print(plot_one(names(res)[i], res[[i]]))
dev.off()
