library(dplyr)
b = import('base')

read_one = function(fname) {
    message(fname)

    extract_flanking = function(...) {
        cols = unlist(list(...))
        cols = na.omit(cols[grepl("Flanking Gene", names(cols))])
        data_frame(gene = b$grep("(^[^ ]+)", cols),
                   ensembl_gene_id = b$grep("(ENS[A-Z0-9]+)", cols),
                   known_cancer = grepl("Cancer Yes", cols))
    }

    cis = readxl::read_xls(fname)[-1,] %>%
        transmute(sample = gsub("[:alnum:_]", "", `Tumour ID`),
                  chr = Chromosome,
                  strand = `Hit Strand`,
                  gene_name = `Hit Ensembl Gene`,
                  ensembl_gene_id = X__1 %catch% NA,
                  known_cancer = grepl("Cancer 1", X__2) %catch% NA,
                  flanking = purrr::pmap(., extract_flanking),
                  assembly = `Assembly Version`,
                  reads = `Read Coverage`)
}

files = list.files("cis_per_tumor", full.names=TRUE)
cis = lapply(files, read_one) %>%
    dplyr::bind_rows()
save(cis, file="cis_per_tumor.RData")
