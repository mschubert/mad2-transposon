library(dplyr)

index_run = function(fastq_dir) {
    regex = "([0-9A-Z]+_[0-9]+)_(S[0-9]+)_R(1)_[0-9]+\\.fastq\\.gz"
    fnames = list.files(fastq_dir, regex, full.names=TRUE)
    keys = sub(regex, "\\1", basename(fnames))
    Rs = length(unique(sub(regex, "\\3", basename(fnames))))
    ends = c("1"="single", "2"="paired")[as.character(Rs)]
#    if (length(ends) != 1)
#        stop("single and paired ends in ", fastq_dir)

    df = tibble(keys, fnames) %>%
        group_by(keys) %>%
        summarize(fastq = list(fnames)) %>%
        mutate(counts = file.path(sub("^seq_fastq", "seq_aligned", fastq_dir),
                                  paste0(keys, ".ReadsPerGene.out.tab"))) %>%
        rowwise() %>%
        mutate(both = list(list(fastq=fastq, counts=counts)))

    setNames(df$both, df$keys)
}

runs = list.files("seq_fastq", "^[0-9A-Z_]+$", all.files=TRUE, full.names=TRUE)
samples = lapply(runs, index_run) %>% setNames(basename(runs))

yaml::write_yaml(samples, file="fastq_index.yaml")
