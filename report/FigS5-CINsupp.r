library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpmisc)
sys = import('sys')
tcga = import('data/tcga')

isCINsig_plot = function(dset) {
    cins = tidyr::gather(dset, "measure", "value", aneup_log2seg, CIN70_Carter2006, HET70)
    stats = cins %>% group_by(measure) %>%
        summarize(res = list(lm(value ~ rev48_stat1_over_wt) %>% broom::tidy() %>%
                             filter(term == "rev48_stat1_over_wt"))) %>%
        tidyr::unnest_wider(res)

    ggplot(cins, aes(x=rev48_stat1_over_wt, y=value)) +
        geom_point(aes(fill=type, shape=factor(p53_mut)), size=2, alpha=0.6) +
        geom_smooth(method="lm", se=FALSE) +
        facet_wrap(~measure, scales="free") +
        geom_label_npc(data=stats, aes(label=sprintf("p=%.2g", p.value)),
                       npcx=0.08, npcy=0.95, size=3, label.size=NA, fill="#ffffffc0") +
        scale_shape_manual(values=c("0"=21, "1"=23), name="TP53 mutation") +
        scale_fill_discrete(name="PAM50 subtype", na.value="#ffffff00") +
        guides(fill = guide_legend(override.aes = list(shape=21))) +
        xlab("CIN + STAT1 ko signature (rev48_stat1_over_wt)")
}

sgl_plot = function(dset) {
    sgl_one = function(x, y) {
        fml = as.formula(sprintf("`%s` ~ `%s`", y, x))
        list(all = lm(fml, data=dset) %>% broom::tidy(),
             wt = lm(fml, data=dset %>% filter(p53_mut == 0)) %>% broom::tidy(),
             mut = lm(fml, data=dset %>% filter(p53_mut == 1)) %>% broom::tidy()) %>%
        bind_rows(.id="p53_status") %>%
        filter(grepl(x, term, fixed=TRUE))
    }
    sgl = expand.grid(x = c("type", "p53_mut", "myc_copy", "MYC", "rev48_stat1_over_wt"),
                      y = c("MYC", "Myc Targets V1"), stringsAsFactors=FALSE) %>%
        filter(x != y) %>%
        rowwise() %>%
            mutate(res = list(sgl_one(x, y))) %>%
        ungroup() %>%
        tidyr::unnest(res) %>%
        filter(!is.na(estimate)) %>%
        mutate(p53_status = factor(p53_status, levels=c("all", "wt", "mut")) %>%
                    forcats::fct_recode("all samples"="all", "p53 wt"="wt", "p53 mut"="mut"))

    ggplot(sgl, aes(x=term, y=-log10(p.value))) +
        geom_col(aes(fill=c("1"="positive","-1"="negative")[as.character(sign(estimate))])) +
        facet_grid(y ~ p53_status, scales="free_x", space="free_x") +
        scale_fill_manual(values=c("positive"="#d6604d", "negative"="#92c5de")) +
        theme(axis.text.x = element_text(hjust=1, angle=45)) +
        labs(title="Single associations of Predictor with Myc/Myc targets", x="Predictor", fill="Regression slope")
}

aov_plot = function(dset) {
    aov_one = function(fml) {
        do_aov = function(mod) car::Anova(mod) %>% broom::tidy() %>%
            mutate(frac_sq = sumsq / sum(rev(sumsq)[-1], na.rm=TRUE),
                   frac_sq_total = sumsq / sum(sumsq, na.rm=TRUE),
                   total_sq = sum(sumsq, na.rm=TRUE) / nobs(mod))
        list(all = lm(fml, dset) %>% do_aov(),
             wt = lm(fml, dset %>% filter(p53_mut == 0)) %>% do_aov(),
             mut = lm(fml, dset %>% filter(p53_mut == 1)) %>% do_aov()) %>%
        bind_rows(.id = "p53_status") %>%
            mutate(rel_sq = total_sq / max(total_sq, na.rm=TRUE))
    }
    anov = list(MYC = aov_one(MYC ~ type + p53_mut + myc_copy + rev48_stat1_over_wt),
                MycV1 = aov_one(`Myc Targets V1` ~ type + p53_mut + MYC + myc_copy + rev48_stat1_over_wt)) %>%
        bind_rows(.id="y") %>%
        mutate(p53_status = factor(p53_status, levels=c("all", "wt", "mut")) %>%
                    forcats::fct_recode("all samples"="all", "p53 wt"="wt", "p53 mut"="mut"),
               y = factor(y, levels=c("MYC", "MycV1", "aneup_log2seg")))

    ggplot(anov, aes(x=rel_sq/2, fill=term, y=frac_sq_total, width=rel_sq)) +
        geom_col() +
        coord_polar("y", start=0) +
        facet_grid(y ~ p53_status) +
        scale_fill_manual(values=c(MYC="#e31a1c", myc_copy="#fb9a99", p53_mut="#35978f",
                                   rev48_stat1_over_wt="#6a3d9a", type="#b15928", Residuals="#dedede")) +
        theme(axis.title = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.grid = element_blank(),
              panel.background = element_blank()) +
        labs(title="Variance explained (Type II ANOVA)", fill="Predictor")
}

purplot = function(dset) {
    ds = dset %>%
        left_join(tcga$purity_estimate() %>% select(Sample, stromal_score, immune_score)) %>%
        mutate(stromal_score = scale(stromal_score),
               immune_score = scale(immune_score),
               purity = scale(purity),
               p53_mut = as.character(p53_mut),
               aneup = scale(aneup_log2seg),
               MYCg = cut(MYC, c(-Inf,median(MYC),Inf)),
               MYC = scale(MYC),
               MycV1 = cut(`Myc Targets V1`, c(-Inf,0,Inf), labels=FALSE),
               MycV1 = forcats::fct_recode(factor(MycV1), "Myc Targets low"="1", "Myc Targets high"="2"))

    ds2 = ds %>%
        tidyr::gather("type", "score", purity, immune_score, stromal_score) %>%
        mutate(type = factor(sub("_score", "", type),
                             levels=c("purity", "immune", "stromal")))

    tests = ds2 %>%
        group_by(type, MycV1, iclass) %>%
            summarize(res = list(broom::tidy(lm(score ~ p53_mut)))) %>%
        ungroup() %>%
        tidyr::unnest(res) %>%
        filter(term == "p53_mut1") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"),
               signif = ifelse(adj.p < 0.1, "FDR<0.1", "n.s."),
               y = ifelse(grepl("low", MycV1), 3, -3))

    mloCIN = ds2 %>% filter(MycV1 == "Myc Targets low", type == "purity", iclass != "noCIN")
    mhiCIN = ds2 %>% filter(MycV1 == "Myc Targets high", type == "purity", iclass != "noCIN")
    imm = ds2 %>% filter(type == "immune", iclass == "CIN_stat1ko")
    t2 = list(
        low_wt = lm(score ~ iclass, data=mloCIN %>% filter(p53_mut == 0)),
        low_mut = lm(score ~ iclass, data=mloCIN %>% filter(p53_mut == 1)),
        high_both = lm(score ~ iclass, data=mhiCIN),
        cmp_wt = lm(score ~ MycV1, data=imm %>% filter(p53_mut == 0)),
        cmp_mut = lm(score ~ MycV1, data=imm %>% filter(p53_mut == 1))
    ) %>%
        lapply(broom::tidy) %>% bind_rows(.id="cmp") %>% filter(term != "(Intercept)") %>%
        mutate(adj.p = p.adjust(p.value, method="fdr"),
               type = factor(rep("immune", 5), levels=levels(ds2$type)), # x
               score = c(-3.8, -2.8, 3, -4, -3.5), # y
               MycV1 = factor(sprintf("Myc Targets %s", c("low", "low", "high", "low", "low")), levels=levels(ds2$MycV1)),
               iclass = factor(c("CIN", "CIN", "CIN", "CIN_stat1ko", "CIN_stat1ko"), levels=levels(ds2$iclass)),
               signif = ifelse(adj.p < 0.1, "FDR<0.1", "n.s."),
               hj = c(0.5, 0.5, 0.5, 1.5, -0.5))

    ggplot(ds2, aes(x=type, fill=type, color=p53_mut, y=score)) +
        geom_boxplot(outlier.shape=NA, position="dodge") +
        facet_grid(MycV1 ~ iclass) +
        scale_color_manual(values=c("0"="#ababab", "1"="#131313")) +
        geom_text(data=tests, aes(y=y, label=sprintf("p=%.2g", p.value), alpha=signif),
                  color="black", size=3, angle=20) +
        geom_text(data=t2, aes(label=sprintf("p=%.2g", p.value), alpha=signif, hjust=hj),
                  color="black", size=3) +
        scale_alpha_manual(values=c("FDR<0.1"=1, "n.s."=0.6)) +
        labs(title = "BRCA by survival class, Myc Targets and p53 status",
             x = "Sample composition or Myc gene expression",
             y = "z-score",
             fill = "Measure",
             alpha = "Difference",
             color = "p53 mut") +
        theme(axis.text.x = element_text(angle=30, hjust=1))
}

sys$run({
    args = sys$cmd$parse(
        opt('b', 'brca', 'rds', '../tcga_myc/dset.rds'),
        opt('p', 'plotfile', 'pdf', 'FigS5-CINsupp.pdf')
    )

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

    p1 = isCINsig_plot(dset)
    p2 = purplot(dset)
    right = sgl_plot(dset) + aov_plot(dset) + plot_layout(nrow=2, heights=c(2,3))
    left = p1 + p2 + plot_layout(heights=c(1,2), ncol=1)
    p = (left | right) + plot_layout(widths=c(1.6,1)) +
        plot_annotation(tag_levels='a') &
        theme(panel.background = element_rect(fill="#f5f5f5"),
              plot.tag = element_text(size=18, face="bold"))

    pdf(args$plotfile, 16, 9)
    print(p)
    dev.off()
})
