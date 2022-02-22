library(dplyr)
library(tidygraph)
sys = import('sys')
idmap = import('process/idmap')
tcga = import('data/tcga')

stat1int_mut = function() {
    op = OmnipathR::import_all_interactions()
    net = OmnipathR::interaction_graph(op) %>%
        as_tbl_graph() %>%
            convert(to_undirected, .clean=TRUE) %>%
            convert(to_simple, .clean=TRUE) %E>%
        select(-.orig_data)
    ints = igraph::neighbors(net, "STAT1")$name
    tot = igraph::V(net)$name

    load_fl = function(coh) tcga$mutations(coh) %>%
        filter(Variant_Classification != "Silent") %>%
        transmute(cohort=coh, Sample=Sample, gene=Hugo_Symbol) %>%
        distinct()
    m2 = lapply(tcga$cohorts(), load_fl) %>%
        dplyr::bind_rows()
    m = m2 %>%
        group_by(cohort, Sample) %>%
            summarize(PIAS1 = length(intersect(gene, "PIAS1")),
                      HEG1 = length(intersect(gene, "HEG1")),
                      SCAF11 = length(intersect(gene, "SCAF11")),
                      EP300 = length(intersect(gene, "EP300")),
                      LPP = length(intersect(gene, "LPP")),
                      TRAPPC9 = length(intersect(gene, "TRAPPC9")),
                      GLIS3 = length(intersect(gene, "GLIS3")),
                      MBNL1 = length(intersect(gene, "MBNL1")),
                      TP53 = length(intersect(gene, "TP53")),
                      EGFR = length(intersect(gene, "EGFR")),
                      TTN = length(intersect(gene, "TTN")),
                      stat1 = length(intersect(gene, "STAT1")),
                      ints = length(intersect(gene, ints)),
                      tot = length(intersect(gene, tot))) %>%
        ungroup() %>%
        mutate(frac = ints/tot) %>%
        arrange(-frac)

    genes = c("Stat1", "Pias1", "Nfkb1", "Ptpn2", "Heg1", "Scaf11", "Ep300", "Trappc9", "Glis3", "Mbnl1", "Trp53", "Myc")
    gex = readRDS("../data/rnaseq/assemble.rds")$expr
    rownames(gex) = idmap$gene(rownames(gex), to="external_gene_name")
    gex = gex[rownames(gex) %in% genes,] %>%
       reshape2::melt() %>%
       as_tibble() %>%
       dplyr::rename(gene=Var1, sample=Var2, expr=value) %>%
       inner_join(meta)
    ggplot(gex, aes(x=gene, group=sample, y=expr)) +
        ggbeeswarm::geom_quasirandom(aes(color=type))
    gdsc = import('data/gdsc')
    ge = gdsc$basal_expression()
    ge = ge[rownames(ge) %in% toupper(genes),] %>%
       reshape2::melt() %>%
       as_tibble() %>%
       dplyr::rename(gene=Var1, sample=Var2, expr=value) %>%
       mutate(cline = gdsc$cosmic$id2name(sample),
              tissue = gdsc$cosmic$id2tissue(sample),
              is_BRCA = tissue == "BRCA",
              label = ifelse(cline %in% c("BT-549", "MCF7", "MDA-MB-231", "A549"), cline, NA))
    ggplot(ge, aes(x=gene, y=expr)) +
        ggbeeswarm::geom_quasirandom(aes(alpha=is_BRCA, color=label, size=!is.na(label))) +
        scale_alpha_manual(values=c("TRUE"=0.8, "FALSE"=0.02)) +
        scale_color_brewer(palette="Set1", na.value="black") +
        scale_size_manual(values=c("TRUE"=3, "FALSE"=1)) #+
#        ggrepel::geom_label_repel(aes(label=label), size=2)

    xx = tcga$aneuploidy() %>%
        select(Sample, aneup=aneup_log2seg) %>%
        inner_join(tcga$purity() %>% select(Sample, cohort, purity=estimate)) %>%
        filter(substr(Sample, 14, 16) == "01A",
               !is.na(aneup) & !is.na(purity)) %>%
        mutate(aneup = aneup / purity,
               aneup_class = cut(aneup, c(0,0.2,Inf))) %>%
        inner_join(m)

    # cands: LPP, GLIS3
    ggplot(xx %>% filter(TP53==1, tot >= 10), aes(x=aneup_class, y=1/tot, color=aneup_class)) +
        ggbeeswarm::geom_quasirandom(alpha=0.5, size=5) +
        scale_color_manual(values=setNames(c("#cc9933", "#5500aa"), c("(0,0.2]", "(0.2,Inf]")), name="Aneuploidy")

    ggplot(xx, aes(x=TRAPPC9, y=tot, color=aneup_class, fill=aneup_class)) +
        geom_point(alpha=0.2) +
        geom_smooth(method="lm") +#, formula=y~s(x, k=10)) +
        scale_x_continuous(trans="log1p", breaks=c(0,1)) +
        scale_y_continuous(trans="log1p", breaks=c(1,10,50,200,500,4000)) +
        scale_color_manual(values=setNames(c("#cc9933", "#5500aa"), c("(0,0.2]", "(0.2,Inf]")), name="Aneuploidy") +
        scale_fill_manual(values=setNames(c("#cc9933", "#5500aa"), c("(0,0.2]", "(0.2,Inf]")), name="Aneuploidy") +
        labs(y="Total number of mutations")

    yy = xx %>%
        filter(stat1 == 0, ints %in% 1:3) %>%
        left_join(m2 %>% filter(gene %in% ints)) %>%
        group_by(gene) %>%
            summarize(n_eup = n_distinct(Sample[aneup_class == "(0,0.1]"]),
                      n_aneup = n_distinct(Sample[aneup_class == "(0.1,Inf]"])) %>%
        arrange(-n_eup)

    zz = inner_join(mouse_ints %>% mutate(gene = toupper(external_gene_name)), yy) %>%
        filter(size >= 3, p.value < 0.1, n_eup+1 >= n_aneup) %>%
        arrange(p.value)

}

sys$run({
    args = sys$cmd$parse(
        opt('a', 'aset', 'rds', '../ploidy_compare/analysis_set.rds'),
        opt('e', 'ext', 'rds', '../cis_analysis/ext_gene.rds'),
        opt('b', 'bionet', 'rds', '../cis_analysis/bionet_omnipath.rds'),
        opt('p', 'plotfile', 'pdf', 'stat1-int-muts.pdf')
    )

    meta = readRDS(args$aset)$meta

    # load mouse ins data from supp fig
    bn = readRDS(args$bionet)
    bn$cis_net = bn$cis_net %>% mutate(hub = centrality_hub())
    common = intersect(igraph::V(bn$cis_net)$name,
                       igraph::V(bn$ext_nets$aneuploidy)$name)
    ext = readRDS(args$ext)$aneuploidy %>%
        mutate(circle = external_gene_name %in% common)
    net_with_stats = bn$ext_nets$aneuploidy %>%
        left_join(ext %>% select(name=external_gene_name, p.value, statistic))

    stat1 = igraph::V(net_with_stats)["Stat1"]
    mouse_ints = net_with_stats %N>%
        as.data.frame() %>% as_tibble() %>%
        filter(external_gene_name %in% igraph::neighbors(net_with_stats, stat1)$name) %>%
        arrange(p.value)

    # load tcga data from supp fig



}
