library(dplyr)
b = import('base')

read_one = function(fname) {
    message(fname)

    extract_flanking = function(...) {
        cols = unlist(list(...))
        flank = cols[grepl("Flanking Gene(__[0-9])?$", names(cols))]
        nona = flank[!is.na(flank)]
        data_frame(gene_name = b$grep("(^[^ ]+)", nona),
                   ensembl_gene_id = b$grep("(ENS[A-Z0-9]+)", nona),
                   known_cancer = grepl("Cancer Yes", nona),
                   hit_dist = c(rev(seq_along(na.omit(flank[1:5]))),
                                seq_along(na.omit(flank[6:10]))))
    }

    cis = readxl::read_xls(fname)[-2,] %>%
        filter(`Tumour ID` != "tumourid") %>%
        transmute(sample = `Tumour ID` %>%
                      sub("PB[0-9]+", "", .) %>%
                      gsub("[^A-Za-z0-9]", "", .) %>%
                      tolower() %>%
                      sub("spl", "s", .) %>%
                      sub("thy", "t", .),
                  chr = Chromosome,
                  position = as.integer(`Transposon Integration Site`),
                  strand = `Transposon Ori`,
                  gene_name = b$grep("(^[^ ]+)", `Hit Ensembl Gene`),
                  ensembl_gene_id = X__1 %catch% b$grep("(ENS[A-Z0-9]+)", `Hit Ensembl Gene`),
                  known_cancer = grepl("Cancer 1", X__2) %catch% grepl("Cancer 1", `Hit Ensembl Gene`),
                  flanking = purrr::pmap(., extract_flanking),
                  assembly = `Assembly Version`,
                  reads = ifelse(`Transposon End` == "3P",
                                 as.integer(`Coverage 3'-end`),
                                 as.integer(`Coverage 5'-end`))) %>%
        filter(sample != "tumourid")
}

files1 = list.files("cis_per_tumor", full.names=TRUE)
files2 = list.files("TraDIS_IS_180622/COMBI_CHROMOSOME_mincov20_xls", full.names=TRUE)
files = c(files1, files2)
nested = lapply(files, read_one) %>%
    dplyr::bind_rows() %>%
    group_by(sample) %>%
    mutate(n_ins_smp = n()) %>%
    ungroup()

hits = nested %>%
    select(-flanking) %>%
    mutate(hit_dist = 0)

flanking = nested %>%
    select(-(gene_name:known_cancer)) %>%
    tidyr::unnest()

cis = dplyr::bind_rows(hits, flanking) %>%
    select(sample:gene_name, hit_dist, everything()) %>%
    arrange(sample, -reads)

save(cis, file="cis_per_tumor.RData")
