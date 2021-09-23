library(dplyr)
library(ggplot2)
library(patchwork)
sys = import('sys')
plt = import('plot')

plot_one = function(title, comp) {
    plots = mapply(function(c, n) plt$volcano(c, label_top=30, repel=TRUE) + ggtitle(n),
        c=comp, n=names(comp), SIMPLIFY=FALSE)
    plt$text(title, size=6) / wrap_plots(plots, design="ABD\nECC") +
        plot_layout(heights=c(1,40))
}

args = sys$cmd$parse(
    opt('d', 'diff_expr', 'rds', 'diff_expr.rds'),
    opt('p', 'plotfile', 'pdf', 'diff_expr.pdf')
)

res = readRDS(args$diff_expr)

pdf(16, 14, file=args$plotfile)
for (i in seq_along(res))
    print(plot_one(names(res)[i], res[[i]]))
dev.off()
