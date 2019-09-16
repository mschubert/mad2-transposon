library(dplyr)
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'infile', 'xlsx with results', 'PR-080-TES_RESULTS REPORT.xlsx'),
    opt('o', 'outfile', 'rds', 'dset.rds'),
    opt('p', 'plotfile', 'pdf', 'dset.pdf'))

pdata = readxl::read_xlsx(args$infile, "RESULTS-PLATE ID PL1004", range="B17:S91")
pmat = data.matrix(pdata[-1])
rownames(pmat) = pdata$SAMPLE
phospho = pdata %>%
    tidyr::gather("protein", "value", -SAMPLE)

cdata = readxl::read_xlsx(args$infile, "RESULTS-PLATE ID PL1005", range="A29:U119") %>%
    filter(!grepl("Standard", SAMPLE))
cmat = data.matrix(cdata[-1])
rownames(cmat) = cdata$SAMPLE
cytokine = cdata %>%
    tidyr::gather("protein", "value", -SAMPLE) %>%
    group_by(protein) %>%
        mutate(value = value - mean(value[SAMPLE == "Background"])) %>%
    ungroup() %>%
    filter(SAMPLE != "Background")
