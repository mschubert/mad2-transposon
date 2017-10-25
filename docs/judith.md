# Transposon insertion methods

Hi Michael,

I have attached a thesis of a student from George Vassiliou lab at Sanger, at
page 52 the cis analysis is described. Also included an excel as an example how
the processed data looks like. And.. the RNAseq sample overview :)

Below I copied some methods from other papers (there are two types of
transposons commonly used, sleeping beauty and piggybac, so the method might be
different between these 2). Let me know if need more info.

Cheers, Judith

## Rad et al., 2014 https://www.nature.com/articles/ng.3164#methods

Mapping of insertion sequences to the mouse genome and identification of common
integration sites.  Mapping of integrations to the mouse genome was performed
using SSAHA2. Query sequences were filtered to contain splinkerette primer
sequences that were located in the transposon inverted terminal repeats.
Redundant sequences that arose from the same tumor and mapped to the same
genomic location were 'collapsed' to a unique integration. To identify regions
in the genome that are more frequently hit by transposons than would be
expected by chance (common insertion sites; CISs), nonredundant insertions were
finally subjected to statistical analysis using TAPDANCE analysis27 or a
framework based on Gaussian Kernel Convolution (GKC) as described earlier28.
GKC analysis was performed at different kernel sizes (at each 10-kb window from
30 to 100 kb). Probability density distribution of pooled piggyBactransposon
insertion sites across the mouse genome was carried out using normal kernels
with a bandwidth of 500,000 bp.

## March et al., 2011 

Identification of CISs. CISs were identified using two complementary approaches
described below.  Gaussian kernel convolution. CIMPLR is a GKC method developed
in a previ- ous study to identify CISs in retroviral insertional mutagenesis
screens18. An enhanced version of this method was developed for Sleeping Beauty
screens, which takes into consideration the local density of TA nucleotides
within the genome. Kernel widths of 30 kb and 120 kb were selected for CIS
detection. For each CIS, the genomic location of its kernel peak was used as
the reference point to assign gene annotations. A CIS peak was associated with
a gene either when it lay within the coding region of a gene or was within 150
kb up- or downstream of its nearest gene. If the distance to the nearest gene
was greater than 150 kb, no gene name was assigned. Using 30-kb and 120-kb
kernels, 919 and 641 CISs were reported, respectively (Supplementary Table
2a,b). To sum- marize, CIS predictions from the two independent kernel widths,
a method of combining predictions across the kernels was devised, such that
where CISs from the two kernels overlapped, those CISs with the smallest
footprints on the genome were reported. Using this procedure, a total of 997
cross-scale CISs were predicted (Supplementary Table 2c). Depending on the size
and location of CISs, genes may be associated with multiple CISs. Consequently,
the 997 cross-scales CISs are associated with 867 unique genes.

## Takeda et al., 2015 https://www.nature.com/articles/ng.3175

Identification of common insertion sites.  CISs were identified as previously
reported53. All informative sequence reads were mapped to the B6 mouse genome
(mm9). To detect CISs, a GKC12 method was employed using 15,000, 30,000,
50,000, 75,000, 120,000 and 240,000 kernel widths53. When CISs were detected
over several kernel widths, the CISs were merged and the smallest window size
is reported. All insertional patterns were checked using the IGB6.4.1 (BioViz)
genome browser, and CISs were removed from the list when the insertions were
not located within a gene. For CIS genes that had two CISs located within them,
the lowest P value is reported. For CISs containing two genes, both genes or
only one of the genes (when there was a strong indication that it was the
correct gene) was reported in the final gene lists (Supplementary Table 2). CIS
genes were then annotated using the Mouse Genome Informatics (MGI) database,
and three CIS genes, AC147162.1, RP24-118B12.1 and snoU13, that were not
annotated in MGI were removed for subsequent analysis. This left 834 CIS genes
from the 3 cohorts.

The likelihood of local hopping of the SB transposon is increased within the
chromosome where the transposon concatamer is located56. The donor site of the
T2/Onc2 transposon used in this study was located on chromosome 1. Therefore,
all insertions on chromosome 1 were removed from the data set for subsequent
analysis. Sfi1 is listed as a single-copy gene in the reference mouse genome.
However, it is estimated that the mouse genome actually has 20–30 copies of
Sfi1 (ref. 61). Therefore, it is very likely that insertions in different Sfi1
loci were annotated to a single Sfi1 gene located on chromosome 11 as in the
reference mouse genome. The Sfi1 gene was therefore removed from our lists of
CIS genes.
