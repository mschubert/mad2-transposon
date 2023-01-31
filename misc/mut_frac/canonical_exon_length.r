library(dplyr)

ens106 = AnnotationHub::AnnotationHub()[["AH100643"]]
is_prot = GeneBiotypeFilter("protein_coding")
is_chr = SeqNameFilter(c(1:22,'X','Y'))
tx = transcripts(ens106, filter=c(is_prot, is_chr)) %>%
    plyranges::filter(tx_is_canonical == 1)
ex = exonsBy(ens106)[names(tx)]
names(ex) = sub("-[0-9]+$", "", tx$tx_external_name)
ex = ex[!is.na(names(ex))] %>%
    lapply(. %>% summarize(glen = sum(end(.) - start(.))) %>% as.data.frame()) %>%
    bind_rows(.id="gene")

saveRDS(ex, file="canonical_exon_length.rds")
