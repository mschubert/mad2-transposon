library(dplyr)
library(plyranges)
library(Gviz)
options(ucscChromosomeNames=FALSE)
b = import('base')
io = import('io')
seq = import('seq')
sys = import('sys')

plot_alignment = function() {
}

de_ctgs = io$read_table("../data/rnaseq_imfusion/merged_ctgs.txt", header=TRUE)
ins = io$read_table("../data/rnaseq_imfusion/insertions.txt", header=TRUE)

sjs = list.files("../data/rnaseq_imfusion", "SJ.out.tab", recursive=TRUE, full.names=TRUE)
sjs = sjs[!grepl("_STARpass1", sjs)]
names(sjs) = b$grep("/([0-9]{3}[st])/", sjs)

bams = list.files("../data/rnaseq_imfusion", "alignment.bam", recursive=TRUE, full.names=TRUE)
names(bams) = b$grep("/([0-9]{3}[st])/", bams)

sample = "442s"
insertion = "INS_4"

ins_id = paste(sample, insertion, sep=".")
ci = filter(ins, id==ins_id)
# can also import GTF via rtracklayer
txdb = GenomicFeatures::makeTxDbFromGFF("../data/genome/Mus_musculus.GRCm38.92.gtf", format="gtf")
grtrack = GeneRegionTrack(txdb, name=ci$gene_name)
gtrack = GenomeAxisTrack(name="GRCm38.92")
altrack = AlignmentsTrack(bams[sample], type=c("coverage","sashimi"), name="RNA-seq")
region = filter(genes, external_gene_name == "Erg")

at = AnnotationTrack(start=ci$position, width=10000, chromosome=ci$seqname,
                     strand=ci$strand, genome="GRCm38", name="PB")

plotTracks(list(gtrack, grtrack, at, altrack),
           from=start(region), to=end(region), chromosome=seqnames(region),
           extend.right=10000, extend.left=10000,
           sizes=c(1,2,1,4), cex.main=1,
           main=with(ci, sprintf("%s: %s @ %s", id, feature_type, gene_name)))


