library(dplyr)
library(plyranges)
library(Gviz)
options(ucscChromosomeNames=FALSE)
b = import('base')
io = import('io')
seq = import('seq')
sys = import('sys')

args = sys$cmd$parse(
    opt('i', 'dir', 'imfusion directory', '../data/rnaseq_imfusion'),
    opt('c', 'ctgs', 'imfusion merged CTGs', '../data/rnaseq_imfusion/merged_ctgs.txt'),
    opt('r', 'rna_ins', 'imfusion insertions', '../data/rnaseq_imfusion/insertions.txt'),
    opt('d', 'dna_ins', '', '../cis_analysis/analysis_set.RData'),
    opt('a', 'gtf', 'assembly GTF', '../data/genome/Mus_musculus.GRCm38.92.gtf'),
    opt('g', 'gene', 'gene name to plot', 'Erg'),
    opt('f', 'flank', 'plot flank around gene', '10000'),
    opt('p', 'plotfile', 'PDF to save to', 'Erg.pdf'))

args$flank = as.integer(args$flank)
de_ctgs = io$read_table(args$ctgs, header=TRUE)
rna_ins = io$read_table(args$rna_ins, header=TRUE)

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
region = genes[genes$external_gene_name == args$gene] %>%
    anchor_start() %>% stretch(as.integer(args$flank)) %>%
    anchor_end() %>% stretch(as.integer(args$flank))
aw = width(region) / 20

dna_ins = io$load(args$dna_ins) %>%
    GenomicRanges::makeGRangesFromDataFrame(start.field="position",
        end.field="position", keep.extra.columns=TRUE) %>%
    join_overlap_intersect(region, .) %>%
    as.data.frame() %>%
    mutate(strand = ifelse(strand == "-", -1L, 1L))

rna_smp = rna_ins %>% filter(gene_name == args$gene) %>% pull(sample)
use_samples = intersect(names(bams), c(rna_smp, dna_ins$sample))

pdf(args$plotfile)
for (sample_id in use_samples) {
    message(sample_id)
    tracks = list(gtrack, grtrack)
    sizes = c(1, 2)

    ci_rna = filter(rna_ins, sample==sample_id & gene_name==args$gene)
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

    ftype = ci_rna$feature_type[1] %or% "???"
    main = sprintf("%s: %s @ %s", sample_id, ftype, args$gene)

    plotTracks(tracks, from=start(region), to=end(region), chromosome=seqnames(region),
               sizes=sizes, cex.main=1, main=main)
}
dev.off()
