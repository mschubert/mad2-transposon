library(dplyr)
library(Gviz)
options(ucscChromosomeNames=FALSE)
b = import('base')
io = import('io')
seq = import('seq')
sys = import('sys')

args = sys$cmd$parse(
    opt('d', 'dir', 'imfusion directory', '../data/rnaseq_imfusion'),
    opt('c', 'ctgs', 'imfusion merged CTGs', '../data/rnaseq_imfusion/merged_ctgs.txt'),
    opt('i', 'insertions', 'imfusion insertions', '../data/rnaseq_imfusion/insertions.txt'),
    opt('a', 'gtf', 'assembly GTF', '../data/genome/Mus_musculus.GRCm38.92.gtf'),
    opt('g', 'gene', 'gene name to plot', 'Erg'),
    opt('f', 'flank', 'plot flank around gene', '10000'),
    opt('p', 'plotfile', 'PDF to save to', 'Erg.pdf'))

args$flank = as.integer(args$flank)
de_ctgs = io$read_table(args$ctgs, header=TRUE)
ins = io$read_table(args$insertions, header=TRUE)

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
region = genes[genes$external_gene_name == args$gene]

sample_with_ins = ins %>%
    filter(gene_name==args$gene) %>%
    pull(sample) %>%
    unique()

pdf(args$plotfile)
for (sample_id in sample_with_ins) {
    altrack = AlignmentsTrack(bams[sample_id], type=c("coverage","sashimi"), name="RNA-seq")
    ci = filter(ins, sample==sample_id & gene_name==args$gene)
    at = AnnotationTrack(start=ci$position, width=10000, chromosome=unique(ci$seqname),
                         strand=ci$strand, genome="GRCm38", name="PB")

    plotTracks(list(gtrack, grtrack, at, altrack),
               from=start(region), to=end(region), chromosome=seqnames(region),
               extend.right=args$flank, extend.left=args$flank,
               sizes=c(1,2,1,4), cex.main=1,
               main=with(ci, sprintf("%s: %s @ %s", id, feature_type, gene_name)))
}
dev.off()
