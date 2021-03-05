library(dplyr)
b = import('base')

read_one = function(fname) {
    message("readxl :: ", fname)

    extract_flanking = function(...) {
        cols = unlist(list(...))
        flank = cols[grepl("^Flanking Gene\\.", names(cols))]
        stopifnot(length(flank) == 10)
        nona = flank[!is.na(flank)]
        data_frame(gene_name = b$grep("(^[^ ]+)", nona),
                   ensembl_gene_id = b$grep("(ENS[A-Z0-9]+)", nona),
                   known_cancer = grepl("Cancer Yes", nona),
                   hit_dist = c(rev(seq_along(na.omit(flank[1:5]))),
                                seq_along(na.omit(flank[6:10]))))
    }

    cis = readxl::read_excel(fname)

    cis = cis %>%
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

files0 = list.files("pilot", "\\.xlsx?$", full.names=TRUE)
files1 = list.files("cis_per_tumor", "\\.xlsx?$", full.names=TRUE)
files2 = list.files("TraDIS_IS_180622/COMBI_CHROMOSOME_mincov20_xls", "\\.xlsx?$", full.names=TRUE)
files = c(files0, files1, files2)
lookup = c( # pilot had different tumor IDs in the xls, standardize them
    "1jspbpb" = "145s",
    "2jspbpb" = "155s",
    "3jspbpb" = "157s",
    "4jspbpb" = "180s",
    "8jspbpb" = "184s",
    "23jspbpb" = "212s"
)
nested = lapply(files, read_one) %>%
    dplyr::bind_rows() %>%
    group_by(sample) %>%
    mutate(sample = ifelse(sample %in% names(lookup), lookup[sample], sample),
           n_ins_smp = n()) %>%
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

stopifnot(all(cis$assembly == "GRCm38"))

saveRDS(cis, file="cis_per_tumor.rds")
