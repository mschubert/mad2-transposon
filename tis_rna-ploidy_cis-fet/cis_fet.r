library(dplyr)
library(cowplot)
io = import('io')
sys = import('sys')

#' FET of aneuploidy low/high vs. gene insertions
#'
#' @param df  Data.frame with 2 rows (one for each condition - e.g. mad2 high
#'  and mad2 low) and the columns 'ins_gene' (number of insertions in a gene in
#'  all samples) and 'ins_total' (number of total insertions across all samples)
#' @return    Results of a Fisher's Exact test for this gene
test_gene = function(df) {
    mat = data.matrix(df[,c('ins_gene','ins_total')])
    fisher.test(mat) %>%
        broom::tidy() %>%
        transmute(condition = df$condition[1],
                  ins_cond = df$ins_gene[1],
                  ins_rest = df$ins_gene[2],
                  log_odds = estimate,
                  p.value = p.value)
}

#' Test all genes for enrichment in insertions
#'
#' @param df  Data.frame with the following fields: 'gene' (identifier),
#'  'ins_gene' (number of insertions in this gene), 'ins_total' (number of
#'  insertions total), and 'condition' (which condition was tested for)
#' @return    Data.frame for result of a Fisher's Exact Test for all genes
test_all = function(df) {
    keep = df %>%
        group_by(gene) %>%
        summarize(ins_gene = sum(ins_gene)) %>%
        filter(ins_gene >= 5) %>%
        pull(gene) %>%
        as.character()
    df = df[df$gene %in% keep,]

    df = df %>%
        group_by(gene) %>%
        tidyr::nest() %>%
        mutate(result = purrr::map(data, test_gene)) %>%
        select(-data) %>%
        tidyr::unnest() %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(p.value)
}

#' Plot the frequency of insertions
#'
#' @param result  Data.frame with results of 'test_all' function
#' @param df      Data.frame, same as input for 'test_all'
#' @param title   Plot title (default: empty string)
#' @return  A ggplot2 object with total number of inserts in groups and highest
#'    scoring genes (numbers and FDR)
ins_plot = function(df, result, title="") {
    ins_all = df %>%
        select(condition, ins_total) %>%
        distinct()

    p1 = ggplot(ins_all, aes(x=condition, y=ins_total, fill=condition)) +
        geom_bar(stat="identity") +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        guides(fill=FALSE)
    p2 = result %>%
        head(10) %>%
        mutate(gene = factor(gene, levels=unique(gene))) %>%
        tidyr::gather("cond", "n_ins", ins_cond, ins_rest)
    p2$condition[p2$cond == "ins_rest"] = sub("high", "low", p2$condition[p2$cond=="ins_cond"])

    if (any(grepl("Intergenic", p2$gene))) {
        p2$n_ins[p2$gene == "Intergenic"] = p2$n_ins[p2$gene == "Intergenic"] / 50
        levels(p2$gene)[levels(p2$gene) == "Intergenic"] = "Intergenic / 50"
    }

    p2 = ggplot(p2, aes(x=gene, y=n_ins, fill=condition)) +
        geom_bar(stat="identity") +
        geom_text(aes(label = sprintf("%.2g", adj.p)), y=15, angle=45) +
        theme(axis.text.x = element_text(angle=45, hjust=1)) +
        ggtitle(title)
    p_mad2 = plot_grid(p1, p2, ncol=2, rel_widths=c(1,3))
}

sys$run({
    # load sample-level aneuploidy w/ cutoffs
    dset = io$load('aneuploidy_mad2.RData')
    cis = io$read_table('../data/cis/171212 CIS RNA seq samples.tsv', header=TRUE) %>%
        mutate(Gene = stringr::str_trim(as.character(Gene))) %>%
        group_by(Sample, Gene) %>%
        summarize(n_ins = n()) %>%
        ungroup() %>%
        narray::construct(n_ins ~ Gene + Sample, fill=0) %>%
        narray::melt()
    colnames(cis) = c("gene", "sample", "n_ins")

    both = inner_join(dset, cis, by="sample") %>%
        select(-aneup, -mad2)

    prep_mad2 = . %>%
        filter(aneup_class != "high") %>%
        group_by(mad2_class, gene) %>%
        summarize(ins_gene = sum(n_ins)) %>%
        ungroup() %>%
        group_by(mad2_class) %>%
        mutate(ins_total = sum(ins_gene),
               condition = paste0("mad2-", mad2_class)) %>%
        ungroup()

    mad2 = prep_mad2(both)
    result_mad2 = test_all(mad2)

#    mad2_coding = prep_mad2(both %>% filter(gene != "Intergenic"))
#    result_mad2_coding = test_all(mad2_coding)

    prep_aneup = . %>%
        filter(mad2_class == "low") %>%
        group_by(aneup_class, gene) %>%
        summarize(ins_gene = sum(n_ins)) %>%
        ungroup() %>%
        group_by(aneup_class) %>%
        mutate(ins_total = sum(ins_gene),
               condition = paste0("aneup-", aneup_class)) %>%
        ungroup()

    aneup = prep_aneup(both)
    result_aneup = test_all(aneup)

#    aneup_coding = prep_aneup(both %>% filter(gene != "Intergenic"))
#    result_aneup_coding = test_all(aneup_coding)

    pdf("cis_fet.pdf")
    print(ins_plot(mad2, result_mad2, "Mad2, all genes"))
#    print(ins_plot(mad2_coding, result_mad2_coding, "Mad2, coding genes"))
    print(ins_plot(aneup, result_aneup, "Aneuploidy, all genes"))
#    print(ins_plot(aneup_coding, result_aneup_coding, "Aneuploidy, coding genes"))
    dev.off()

    save(mad2, aneup, file="cis_fet.RData")
})
