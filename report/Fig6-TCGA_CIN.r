library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(patchwork)
library(ggpmisc)
sys = import('sys')
tcga = import('data/tcga')

#todo: add what is now Fig4 BRCA sig compare

#todo: if instead of stat1ko, put in aneup/CIN70/E2F? @supp
# + naive assocs with those @supp

pancan_myc_stat = function() {
#    pur = tcga$purity() %>% select(Sample, cohort, purity=estimate) #FIXME: why is xcell different?
    pur = tcga$purity_estimate()
    cohorts = unique(pur$cohort)
    cin70 = lapply(cohorts, function(c) tcga$gsva(c, "CIN")["CIN70_Carter2006",]) %>%
        do.call(c, .) %>% stack() %>% as_tibble() %>% select(Sample=ind, CIN70=values)
    myc = lapply(cohorts, function(c) tcga$gsva(c, "MSigDB_Hallmark_2020")["Myc Targets V1",]) %>%
        do.call(c, .) %>% stack() %>% as_tibble() %>% select(Sample=ind, MycV1=values)
    stat = lapply(cohorts, function(c) tcga$gsva(c, "DoRothEA")["STAT1 (a)",]) %>%
        do.call(c, .) %>% stack() %>% as_tibble() %>% select(Sample=ind, STAT1=values)
    p53mut = lapply(cohorts, function(c) tcga$mutations(c) %>% filter(Hugo_Symbol=="TP53") %>%
                    select(Sample) %>% mutate(p53_mut=1)) %>% bind_rows()

    ds = tcga$aneuploidy() %>% select(Sample, aneup=aneup_log2seg) %>%
        filter(substr(Sample, 14, 16) == "01A") %>%
        inner_join(pur) %>%
        inner_join(cin70) %>%
        inner_join(myc) %>%
        inner_join(stat) %>%
        left_join(p53mut) %>%
        mutate(p53_mut = ifelse(is.na(p53_mut), 0, p53_mut),
               aneup = aneup / purity,
               cohort = tcga$barcode2study(Sample))
    dsRep = bind_rows(list(
        ds %>% mutate(cohort2 = cohort),
        ds %>% filter(p53_mut == 0) %>% mutate(cohort2 = paste0(cohort, "_wt")),
        ds %>% filter(p53_mut == 1) %>% mutate(cohort2 = paste0(cohort, "_mut"))
    ))

    x= dsRep %>% group_by(cohort2) %>%
        summarize(res = list(broom::tidy(lm(STAT1 ~ purity + aneup)))) %>%
        tidyr::unnest(res)
    y= dsRep %>% group_by(cohort2) %>%
        summarize(res = list(broom::tidy(lm(MycV1 ~ purity + aneup)))) %>%
        tidyr::unnest(res)

    ds2 = inner_join(
        x %>% filter(term == "aneup") %>% select(cohort2, STAT1=estimate, seSTAT1=std.error),
        y %>% filter(term == "aneup") %>% select(cohort2, MycV1=estimate, seMycV1=std.error)
    ) %>% inner_join(dsRep %>% group_by(cohort, cohort2) %>% summarize(n=n())) %>%
        mutate(meanSE = (seSTAT1 + seMycV1)/2)
    ggplot(ds2 %>% filter(cohort2!="KIRC_mut"), aes(x=STAT1, y=MycV1, color=cohort)) +
        geom_point(aes(size=n, alpha=meanSE)) +
        scale_alpha_continuous(trans="reverse") +
#        scale_color_brewer(palette="Set3") +
        scale_size_area() +
        ggrepel::geom_text_repel(aes(label=cohort2)) +
        geom_vline(xintercept=0, linetype="dashed") +
        geom_hline(yintercept=0, linetype="dashed")

    #todo: could show that Myc on avg increases, stat1 unchanged -> separate BRCA in STAT1+/- via KO sig
    # -> show that acuteCIN-high and KO-high are more aneup then no-CIN + effect in surv
}

stat1int_mut = function() {
    library(tidygraph)
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
        transmute(cohort=coh, Sample=Sample, gene=Hugo_Symbol)
    m = lapply(tcga$cohorts(), load_fl) %>%
        dplyr::bind_rows() %>%
        group_by(cohort, Sample) %>%
        summarize(stat1 = length(intersect(gene, "STAT1")),
                  ints = length(intersect(gene, ints)),
                  tot = length(intersect(gene, tot))) %>%
        mutate(frac = ints/tot) %>%
        arrange(-frac)

    xx= inner_join(m ,ds) %>% filter(!is.na(aneup)) %>% #%>% left_join(dset %>% select(Sample, rev48_stat1_over_wt))
        mutate(aneup_class=cut(aneup, c(0,0.1,Inf)))
    ggplot(xx %>% filter(stat1==1), aes(x=aneup_class, y=1/tot, color=aneup_class)) +
        ggbeeswarm::geom_quasirandom(alpha=0.5, size=5) +
        ggsignif::geom_signif(comparisons=list(c("(0,0.1]", "(0.1,Inf]"))) +
        labs(y="STAT1 as fraction of mutated genes")
    ggplot(xx, aes(x=ints, y=tot, color=aneup_class, fill=aneup_class)) +
        geom_point(alpha=0.2) +
        geom_smooth(method="gam", formula=y~s(x, k=10)) + #todo: can I test smooth difference using mgcv?
#        facet_wrap(~cohort) + # set k=5 above
        scale_x_continuous(trans="log1p") + scale_y_continuous(trans="log1p")
    m1 = mgcv::gam(log10(tot+1) ~ s(log10(ints+1), k=10), data=xx)
    m2 = mgcv::gam(log10(tot+1) ~ aneup_class + s(log10(ints+1), k=10), data=xx)
    broom::tidy(anova(m1, m2, test="F"))

    #todo: find IFNA loss assocs first & link them to aneup/STAT1ko sig
    #todo: add RUBIC BEM incl IFNA loss, then split out IFNa loss + new marker(s)?


    # [x] mut mgcv test diff
    # [ ] stat1ko split survival for all TCGA cohorts
    # [ ] stat interactor copy number
}
stat1int_cna = function() {
    cna_bem = readr::read_tsv("~/Downloads/Tumours_CNV_BEMs/PANCAN_CNA_BEM.rdata.txt")
    # AQUAE_p_TCGA_112_304_b2_N_GenomeWideSNP_6_D06_1348308 etc.
}

mut_stat1ko = function() {
    # which mutations are associated with stat1ko signature score? [bem?]

    bem = tcga$bem()
    yy = ds %>% mutate(Sample = substr(Sample, 1, 12)) %>% filter(cohort == "BRCA")
    narray::intersect(bem, yy$Sample, along=1)
    bem = bem[,colSums(bem) >= 4]

    s1 = yy$STAT1 #todo: actual stat1 ko sig
    purity = yy$purity
    res = st$lm(s1 ~ purity + bem) %>% as_tibble() %>%
        filter(term == "bem") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value) %>% filter(adj.p < 0.2)


    bem = tcga$bem()
    bem = bem[,colSums(bem) >= 4]
    rownames(bem) = paste0(rownames(bem), "-01A")
    narray::intersect(dset$Sample, bem, along=1)
    purity = dset$purity
    s1 = dset$rev48_stat1_over_wt
    aneup = dset$aneup_log2seg
    type = dset$type
    res = st$lm(s1 ~ purity + type + aneup + bem) %>%
        as_tibble() %>%
        filter(term == "bem") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value) #%>% filter(adj.p < 0.2)

    s2 = dset$wt_rev48_over_dmso
    res2 = st$lm(s2 ~ purity + type + aneup + bem) %>%
        as_tibble() %>%
        filter(term == "bem") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr")) %>%
        arrange(adj.p, p.value) #%>% filter(adj.p < 0.2)

    inner_join(res %>% select(bem, size, e1=estimate, p1=p.value, s1=statistic),
               res2 %>% select(bem, e2=estimate, p2=p.value, s2=statistic)) %>%
        mutate(pm = pmin(p1, p2)) %>%
        ggplot(aes(x=s1, y=s2)) +
            geom_point(aes(size=size, alpha=-log10(pm))) +
            ggrepel::geom_text_repel(aes(label=ifelse(s1>2|s2< -3.5|s1-s2>3.5, bem, NA)), size=2) +
            geom_smooth(method="lm") +
            geom_hline(yintercept=0) + geom_vline(xintercept=0) +
            labs(x="stat1ko", y="acute cin")
}

calc_surv = function(dset) {
    calc_one = function(field) {
        fs = rlang::sym(field)
        dset2 = dset %>%
            mutate(field = case_when(
            wt_rev48_over_dmso<0 & !! fs>0 ~ "CIN_field",
            wt_rev48_over_dmso>0 ~ "CIN",
            wt_rev48_over_dmso<0 & !! fs<0 ~ "noCIN"
        ), field=relevel(factor(field), "noCIN"))
        fml = Surv(OS_years, vital_status) ~ age_at_diagnosis + purity + field
        list(#all = coxph(fml, data=dset2),
             wt = coxph(fml, data=dset2 %>% filter(p53_mut == 0)),
             mut = coxph(fml, data=dset2 %>% filter(p53_mut == 1))) %>%
            lapply(broom::tidy) %>% bind_rows(.id="p53_status")
    }

    fields = c("rev48_stat1_over_wt", "CIN70_Carter2006", "HET70",
               "Myc Targets V1", "E2F Targets")
    sapply(fields, calc_one, simplify=FALSE) %>% bind_rows(.id="assoc") %>%
        filter(term == "fieldCIN_field") %>% select(-term)
}

survplot = function(dset) {
    p53wt = dset %>% filter(p53_mut == 0)
    p53mut = dset %>% filter(p53_mut != 0)

    iclass_cmp = paste0("iclass", rev(levels(dset$iclass))[1])
    m1 = coxph(Surv(OS_years, vital_status) ~ age_at_diagnosis + purity + iclass, data=p53wt)
    m1p = broom::tidy(m1) %>% filter(term == iclass_cmp) %>% pull(p.value)
    m2 = coxph(Surv(OS_years, vital_status) ~ age_at_diagnosis + purity + iclass, data=p53mut)
    m2p = broom::tidy(m2) %>% filter(term == iclass_cmp) %>% pull(p.value)

    pal = c("#ababab", "blue", "#ad07e3")
    lab = c("No CIN signature", "Acute CIN signature", "STAT1ko CIN signature")

    fit1 = survfit(Surv(OS_years, vital_status) ~ iclass, data=p53wt)
    ps1 = ggsurvplot(fit1, data=p53wt, xlim=c(0,10), break.time.by=2.5, palette=pal, legend.labs=lab)$plot +
        ylim(c(0.25,1)) +
        labs(x = "Overall survival (years)",
             title = "Acute CIN response vs. STAT1 ko",
             subtitle = sprintf("p53 wt (n=%i)", sum(fit1$n))) +
        annotate("text_npc", npcx=0.1, npcy=0.1,
                 label=sprintf("CIN STAT1ko vs. no CIN p=%.2g\nCIN70 n.s.\nMyc Targets V1 n.s.\nE2F Targets n.s.", m1p))

    fit2 = survfit(Surv(OS_years, vital_status) ~ iclass, data=p53mut)
    ps2 = ggsurvplot(fit2, data=p53mut, xlim=c(0,10), break.time.by=2.5, palette=pal, legend.labs=lab)$plot +
        ylim(c(0.25,1)) +
        labs(x = "Overall survival (years)",
             subtitle = sprintf("p53 mut (n=%i)", sum(fit2$n))) +
        annotate("text_npc", npcx=0.1, npcy=0.1,
                 label=sprintf("CIN STAT1ko vs. no CIN p=%.2g\nCIN70 n.s.\nMyc Targets V1 n.s.\nE2F Targets n.s.", m2p))

#    other = calc_surv(dset)
#    ggplot(other, aes(x=assoc, fill=p53_status, y=-log10(p.value), alpha=p.value<0.05)) +
#        geom_col(position="dodge") +
#        geom_hline(yintercept=-log10(0.05), linetype="dashed") +
#        geom_text(aes(label=sprintf("  %.2g  ", p.value)), position=position_dodge(width=1),
#                  hjust=ifelse(other$p.value<0.3, 1, 0)) +
#        scale_alpha_manual(values=c("TRUE"=0.9, "FALSE"=0.4)) +
#        coord_flip()

    ps1 + ps2 + plot_layout(guides="collect") & theme(legend.direction = "vertical")
}

sys$run({
    args = sys$cmd$parse(
        opt('b', 'brca', 'rds', '../tcga_myc/dset.rds'),
        opt('p', 'plotfile', 'pdf', 'Fig6-TCGA_CIN.pdf')
    )

    scde = readRDS("../data/scRNA_cancer/dset.rds")

    brca = readRDS(args$brca)
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

    survplot(dset)
})

surv_cohort = function(cohort) {
#    clin = tcga$clinical(cohort) %>%
#        transmute(patient = submitter_id,
#                  vital_status = as.integer(factor(vital_status)) - 1,
#                  age_at_diagnosis = age_at_diagnosis / 365,
#                  OS_years = pmax(days_to_death, days_to_last_follow_up, na.rm=TRUE) / 365,
#                  vital_status = ifelse(OS_years > 10, 0, vital_status),
#                  OS_years = pmin(OS_years, 10))
    eset = tcga$rna_seq(cohort, trans="vst", annot=TRUE)
    clin = SummarizedExperiment::colData(eset) %>% as.data.frame() %>% as_tibble() %>%
        transmute(Sample = sample,
                  vital_status = as.integer(factor(vital_status)) - 1,
                  age_at_diagnosis = age_at_diagnosis / 365,
                  OS_years = pmax(days_to_death, days_to_last_follow_up, na.rm=TRUE) / 365,
                  vital_status = ifelse(OS_years > 10, 0, vital_status),
                  OS_years = pmin(OS_years, 10))

    gsv = tcga$gsva(cohort, "CIN")[c("wt_rev48_over_dmso", "rev48_stat1_over_wt"),] %>%
        t() %>% as.data.frame() %>% tibble::rownames_to_column("Sample") %>% as_tibble()
    mut = tcga$mutations(cohort) %>% filter(Hugo_Symbol=="TP53") %>%
        select(Sample) %>% mutate(p53_mut=1)
    pur = tcga$purity() %>% transmute(Sample=Sample, purity=estimate)

    dset2 = gsv %>%
        filter(substr(Sample, 14, 16) == "01A") %>%
        inner_join(clin) %>%
        left_join(mut) %>%
        inner_join(pur) %>%
        mutate(p53_mut = ifelse(is.na(p53_mut), 0, p53_mut))

    dset2 = dset2 %>%
        mutate(iclass = case_when(
            wt_rev48_over_dmso<0 & rev48_stat1_over_wt>0 ~ "CIN_stat1ko",
            wt_rev48_over_dmso>0 ~ "CIN",
            wt_rev48_over_dmso<0 & rev48_stat1_over_wt<0 ~ "noCIN"
        ), iclass=relevel(factor(iclass), "noCIN"))

    survplot(dset)

    x = dset2 %>% inner_join(dset %>% transmute(Sample=Sample, rev2=wt_rev48_over_dmso, stat2=rev48_stat1_over_wt))
    plot(x$wt_rev48_over_dmso, x$rev2) # ok
    plot(x$rev48_stat1_over_wt, x$stat2) # no cor?!
}
