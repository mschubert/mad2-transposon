library(dplyr)
seq = import('seq')
idmap = import('process/idmap')

# filter derived from STAR manual for Chimeric.out.junction
read_file = function(fname) {
    df = read.table(fname, sep="\t") %>%
        filter(V1 != "MT" & V4 != "MT" & V7 > 0 & V8+V9 <= 5) %>%
        group_by(V1, V2, V3, V4, V5, V6, V7, V8, V9) %>%
        summarize(reads = dplyr::n()) %>%
        ungroup() %>%
        arrange(-reads) %>%
        filter(reads >= 2)
}

juncs = list.files(".", pattern="Chimeric.out.junction",
                   recursive=TRUE, full.names=TRUE)
names(juncs) = strsplit(juncs, "/",) %>% lapply(`[`, 2)

coords = seq$coords$gene(dset="mmusculus_gene_ensembl", assembly="GRCm38",
                         granges=TRUE)

support = lapply(juncs, read_file) %>%
    bind_rows(.id="sample") %>%
    mutate(idx = 1:nrow(.))
from = makeGRangesFromDataFrame(support, start.field="V2", end.field="V2",
                                        seqnames.field="V1", strand.field="V3") %>%
    plyranges::mutate(idx = 1:length(.)) %>%
    plyranges::join_overlap_left(coords) %>%
    as.data.frame() %>%
    select(idx, gene1=external_gene_name)
to = makeGRangesFromDataFrame(support, start.field="V5", end.field="V5",
                                      seqnames.field="V4", strand.field="V6") %>%
    plyranges::mutate(idx = 1:length(.)) %>%
    plyranges::join_overlap_left(coords) %>%
    as.data.frame() %>%
    select(idx, gene2=external_gene_name)

merge = support %>%
    left_join(from, by="idx") %>%
    left_join(to, by="idx") %>%
    select(-idx) %>%
    mutate(gene1 = ifelse(V1 == "pA6-GrOnc", "pA6-GrOnc", gene1),
           gene2 = ifelse(V4 == "pA6-GrOnc", "pA6-GrOnc", gene2)) %>%
    filter(!is.na(gene1) & !is.na(gene2))

merge %>% filter(gene1 == "Erg" | gene2=="Erg") %>% arrange(-reads)
