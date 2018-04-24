library(dplyr)
b = import('base')

read_one = function(fname) {
    message(fname)

    extract_flanking = function(...) {
        cols = unlist(list(...))
        cols = na.omit(cols[grepl("Flanking Gene", names(cols))])
        data_frame(gene_name = b$grep("(^[^ ]+)", cols),
                   ensembl_gene_id = b$grep("(ENS[A-Z0-9]+)", cols),
                   known_cancer = grepl("Cancer Yes", cols))
    }

    cis = readxl::read_xls(fname)[-1,] %>%
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
                                 as.integer(`Coverage 5'-end`)))
}

files = list.files("cis_per_tumor", full.names=TRUE)
nested = lapply(files, read_one) %>%
    dplyr::bind_rows() %>%
    group_by(sample) %>%
    mutate(n_ins_smp = n()) %>%
    ungroup()

hits = nested %>%
    select(-flanking) %>%
    mutate(type = "hit")

flanking = nested %>%
    select(-(gene_name:known_cancer)) %>%
    tidyr::unnest() %>%
    mutate(type = "flanking")

cis = dplyr::bind_rows(hits, flanking) %>%
    select(sample:gene_name, type, everything()) %>%
    arrange(sample, -reads)

save(cis, file="cis_per_tumor.RData")
