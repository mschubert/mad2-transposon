library(cowplot)
library(dplyr)
library(plyranges)
library(Gviz)
options(ucscChromosomeNames=FALSE)
b = import('base')
io = import('io')
seq = import('seq')
sys = import('sys')
idmap = import('process/idmap')

args = sys$cmd$parse(
    opt('i', 'dir', 'imfusion directory', '../data/rnaseq_imfusion'),
    opt('e', 'exons', 'expression tsv', '../data/rnaseq_imfusion/exon_counts.txt'),
    opt('c', 'ctgs', 'imfusion merged CTGs', '../data/rnaseq_imfusion/merged_ctgs.txt'),
    opt('r', 'rna_ins', 'imfusion insertions', '../data/rnaseq_imfusion/insertions.txt'),
    opt('d', 'dna_ins', '', '../cis_analysis/analysis_set.RData'),
    opt('m', 'meta', 'metadata', '../ploidy_compare/analysis_set.RData'),
    opt('a', 'gtf', 'assembly GTF', '../data/genome/Mus_musculus.GRCm38.92.gtf'),
    opt('g', 'gene', 'gene name to plot', 'Erg'),
    opt('f', 'flank', 'plot flank around gene', '10000'),
    opt('p', 'plotfile', 'PDF to save to', 'Erg.pdf'))

all_genes = strsplit(args$gene, "+", fixed=TRUE)[[1]]
meta = io$load(args$meta)
de_ctgs = io$read_table(args$ctgs, header=TRUE)
rna_ins = io$read_table(args$rna_ins, header=TRUE)
exons = io$read_table(args$exons, header=TRUE)
expr = DESeq2::DESeqDataSetFromMatrix(exons[,6:ncol(exons)],
        data.frame(id=colnames(exons)[6:ncol(exons)]), ~1) %>%
    DESeq2::estimateSizeFactors() %>%
    DESeq2::counts(normalized=TRUE)
exons[,6:ncol(exons)] = expr
exons = exons %>%
    filter(gene_id %in% idmap$gene(all_genes, from="external_gene_name",
        to="ensembl_gene_id", dset="mmusculus_gene_ensembl")) %>%
    mutate(position = 1:dplyr::n() - 1,
           width = abs(start - end),
           start = (cumsum(width) - width) + position * 100,
           end = cumsum(width) + position * 100) %>%
    select(-gene_id, -chr, -strand) %>%
    tidyr::gather("sample", "reads", -position, -width, -start, -end) %>%
    filter(reads != 0) %>%
    group_by(sample) %>%
    mutate(reads_per_kb = reads / width * 1000,
           next_start = c(start[-1], NA),
           next_reads = c(reads_per_kb[-1], NA),
           label = ifelse(position == max(position), sample, NA)) %>%
    ungroup()

sjs = list.files(args$dir, "SJ.out.tab", recursive=TRUE, full.names=TRUE)
sjs = sjs[!grepl("_STARpass1", sjs)]
names(sjs) = b$grep("/([0-9]{3}[st])/", sjs)

bams = list.files(args$dir, "alignment.bam", recursive=TRUE, full.names=TRUE)
names(bams) = b$grep("/([0-9]{3}[st])/", bams)

# can also import GTF via rtracklayer
txdb = GenomicFeatures::makeTxDbFromGFF(args$gtf, format="gtf")
grtrack = GeneRegionTrack(txdb, name=args$gene)
gtrack = GenomeAxisTrack(name="GRCm38.92")
genes = seq$coords$gene(dset="mmusculus_gene_ensembl", granges=TRUE) # get this out of txdb obj?
region = genes[genes$external_gene_name %in% all_genes] %>%
    anchor_5p() %>% stretch(as.integer(args$flank)) %>%
    anchor_3p() %>% stretch(as.integer(args$flank))
aw = width(region) / 20

dna_ins = io$load(args$dna_ins) %>%
    GenomicRanges::makeGRangesFromDataFrame(start.field="position",
        end.field="position", keep.extra.columns=TRUE) %>%
    join_overlap_intersect(region, .) %>%
    as.data.frame() %>%
    mutate(strand = ifelse(strand == "-", -1L, 1L))

rna_smp = rna_ins %>% filter(gene_name %in% all_genes) %>% pull(sample)
use_samples = c(intersect(names(bams), c(rna_smp, dna_ins$sample)),
                setdiff(names(bams), c(rna_smp, dna_ins$sample)))

exon_ins = exons %>%
    mutate(ins_type = 2*(sample %in% rna_smp) + (sample %in% dna_ins$sample),
           ins_type = factor(ins_type, levels=0:3)) %>%
    inner_join(meta %>% select(sample, type, aneuploidy)) %>%
    mutate(type = ifelse(is.na(type), "unknown", type))
levels(exon_ins$ins_type) = c("none", "DNA", "RNA", "both")

pdf(args$plotfile)
ggplot(exon_ins, aes(color=ins_type, alpha=aneuploidy)) +
    geom_segment(aes(x=end, xend=next_start, y=reads_per_kb, yend=next_reads),
                 size=0.2, linetype="dashed", na.rm=TRUE) +
    geom_segment(aes(x=start, xend=end, y=reads_per_kb, yend=reads_per_kb), size=0.7) +
    scale_color_manual(values=c("#b3b3b3", "#e41a1c", "#377eb8", "#4daf4a"),
                       labels=levels(exon_ins$ins_type)) +
    scale_y_log10() +
    ggrepel::geom_text_repel(aes(x=end+200, y=reads_per_kb, label=label),
                             alpha=1, na.rm=TRUE, size=2, box.padding=0) +
    xlab("distance_from_start") +
    ggtitle(paste(args$gene, "exon expression")) +
    facet_wrap(~type)

for (sample_id in use_samples) {
    message(sample_id)
    tracks = list(gtrack, grtrack)
    sizes = c(1, 2)

    ci_rna = filter(rna_ins, sample==sample_id & gene_name %in% all_genes)
    ci_dna = filter(dna_ins, sample==sample_id)
    chr = unique(c(ci_rna$seqname, ci_dna$seqnames))
    if (nrow(ci_dna) > 0) {
        tracks$dna_pb = AnnotationTrack(width=aw, genome="GRCm38", name="PB DNA",
            chromosome=chr, start=ci_dna$start, strand=ci_dna$strand,
            fill=rep("red", nrow(ci_dna)))
        sizes = c(sizes, 1)
    }
    if (nrow(ci_rna) > 0) {
        tracks$rna_pb = AnnotationTrack(width=aw, genome="GRCm38", name="PB RNA",
            chromosome=chr, start=ci_rna$position, strand=ci_rna$strand,
            fill=rep("blue", nrow(ci_rna)))
        sizes = c(sizes, 1)
    }
    tracks$rna = AlignmentsTrack(bams[sample_id], type=c("coverage","sashimi"), name="RNA-seq")
    sizes = c(sizes, 4)

    ftype = ci_rna$feature_name[1] %or% "???"
    main = sprintf("%s: %s @ %s", sample_id, ftype, args$gene)

    plotTracks(tracks, from=start(region), to=end(region), chromosome=seqnames(region),
               sizes=sizes, cex.main=1, main=main)
}
dev.off()
