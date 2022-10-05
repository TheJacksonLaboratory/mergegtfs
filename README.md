# MERGEGTFS

A script for merging transcript records from multiple 
[GTF2.2](http://mblab.wustl.edu/GTF22.html) formatted input files into 
a single non-redundant GTF2.2 formatted output file. Accommodates various 
local and global rearrangements, including inter-chromosomal fusions. 

---

## SYNTAX

```
usage: mergegtfs.py [-h] [--tol_sj TOL_SJ] [--tol_tss TOL_TSS] [--tol_tts TOL_TTS]
                    [--p_exon_overlap P_EXON_OVERLAP] [--p_exons_overlap P_EXONS_OVERLAP]
                    [--max_intron_length MAX_INTRON_LENGTH] [--sort_exons] [--rev_neg_exons]
                    [--gene_prefix GENE_PREFIX] [--output_prefix OUTPUT_PREFIX]
                    gtf_list_file

Merges redundant transcripts from multiple GTF2.2 formatted files listed in gtf_list_file, resulting in a non-
redundant union gtf.

positional arguments:
  gtf_list_file         Tab-delimited file with 1 unique label followed by 1 GTF2.2 file path per line

optional arguments:
  -h, --help            show this help message and exit
  --tol_sj TOL_SJ       Tolerance (bp) for matching splice junction coordinates (default: 0)
  --tol_tss TOL_TSS     Tolerance (bp) for matching transcript start coordinates (default: 0)
  --tol_tts TOL_TTS     Tolerance (bp) for matching transcript end coordinates (default: 0)
  --p_exon_overlap P_EXON_OVERLAP
                        Minimum proportion overlap between two exons needed for gene matching (default: 0.5)
  --p_exons_overlap P_EXONS_OVERLAP
                        Minimum proportion of overlapping exons needed for gene matching (default: 0.1)
  --max_intron_length MAX_INTRON_LENGTH
                        Maximum expected intron size; used for identifying fusions (default: 1250000)
  --sort_exons          For each transcript, sort exons in ascending order by (begin, end) (default: False)
  --rev_neg_exons       For each transcript, reverse order of negative strand exons (default: False)
  --gene_prefix GENE_PREFIX
                        Prefix for gene_ids and transcript_ids (default: LOC.)
  --output_prefix OUTPUT_PREFIX
                        Prefix for output file names (default: union)

Input GTF2.2 formatted files are required to have exon features with attributes that include 'gene_id' and
'transcript_id'. Output GTF2.2 formatted file includes 'exon' and 'transcript' features, both with attributes
(in order) 'gene_id' and 'transcript_id'. New gene_ids will be assigned in accordance with --p_exon_overlap
(min overlap for matching exons) and --p_exons_overlap (min proportion of matched exons for matching gene_ids).
New transcript_ids are assigned sequentially for each gene. New gene_ids and transcript_ids both have prefix
specified by --gene_prefix. It is assumed that the order of exons in the input GTF file reflects the actual
order in the spliced transcript; it is also assumed that for negative strand transcripts, the last exon is
listed first, and the first exon is listed last. If the order of exon records for negative strand transcripts
in the input GTF file begins with the last exon and ends with the first, make sure to set '--rev_neg_exons'. If
the order of exons in the input file is indeterminate and does not reflect exon order in the transcript, please
set '--sort_exons'; however, this may collapse isoforms where exon order really is shuffled due to e.g. local
genomic rearrangement.
```

---

## INPUTS

Below are the contents of an example `gtf_list_file` specifying 
[GTF2.2](http://mblab.wustl.edu/GTF22.html) formatted files from three samples 
that are to be merged. The format of the file is text with two tab-delimited 
columns and no header row. The first column is a unique sample label, while 
the second column contains a relative or absolute path to the GTF file 
corresponding to that sample. Lines beginning with `#` are treated as 
comments and, along with blank lines, are ignored:

```
# this is a comment; below are sample labels and GTF file paths:
sample1    ~/project/sample1/collapse.gtf
sample2    ~/project/sample2/collapse.gtf
sample3    ~/project/sample3/collapse.gtf
```

#### Gory details

Input GTF2.2 files are expected to have exon records (column 3 is `exon`), 
each of which has 9 tab-delimited columns with attributes (9th column) 
specifying the gene ID and transcript ID associated with the exon in the 
standard format (including the spaces, double-quotes and semicolons):

`gene_id "gene167"; transcript_id "gene167.3";`

Attributes can occur in any order and additional components can be included 
in the attributes as long as the quoting, spaces and semicolons are properly 
used to delimit each component:

`gene_name "Mcl1"; transcript_id "gene167.3"; evidence "3"; gene_id "gene167";`

Accommodating potential genomic rearrangements which give rise to rearranged
transcripts requires certain assumptions about the input GTF file. Exon 
records for each transcript should be listed in the order they are found in 
the final spliced transcript. If you are not interested in genomic 
rearrangements/fusions, and exons for a transcript are not in order, you can 
use the `--sort_exons` flag to disregard exon order. It is also 
assumed that exons on both the positive and negative strands are listed in 
positive strand order. That is, for positive strand transcripts, the first 
exon occuring in the spliced transcript should be listed first, and the last 
exon listed last. By contrast, for negative strand transcripts, the last exon 
will be listed first, and the first exon listed last. One can accommodate GTF 
files with negative strand transcripts exons in negative strand order by 
setting the `--rev_neg_exons` flag. The `--sort_exons` flag will override the 
`--rev_neg_exons` flag if both are specified.

---

## OUTPUTS

Output is composed of a non-redundant GTF2.2 file containing transcript 
records, and a tab-delimited file mapping input identifiers to output 
identifiers. If several transcripts are merged together, one is chosen 
as exemplar and serves as the sole representative of the merged set in 
the output. The mapping file specifies if the transcript appears to be 
rearranged or fused, as well as whether the transcript was used as exemplar 
in the output. 

#### Gory details

By default, the GTF2.2 file is named `union.gtf` and the mapping file is 
named `union.xrefs.tsv`. You can change the prefix `union` using the 
parameter `--output_prefix`.

Transcripts are merged based on the tolerances set with `--tol_tss`, 
`--tol_sj`, and `--tol_tts`. These tolerances are used when comparing 
coordinates of corresponding exons from two transcripts. If the 
strand or the number of exons differs between two transcripts, the transcripts 
will not be merged. If, for any exon, the coordinates of the exon of one 
transcript differ from those of the corresponding exon of the second 
transcript by more than the specified tolerances, those transcripts will not 
be merged. If the strands are the same, the number of exons are the same, and 
all the exon coordinates do not differ by more than the specified tolerances, 
the two transcripts will be merged. When merging, if the start coordinate of 
the first exon of one transcript is before the start coordinate of the first 
exon of a matching transcript, the first transcript will be kept as the 
exemplar. If start coordinates match, the transcript from the first file 
listed in the `gtf_list_file` will serve as exemplar. If both transcripts 
have the same start coordinates and come from the same file, the one which 
occurs first in the file is used as exemplar. Merging is tracked in 
`union.xrefs.tsv`. Whether a transcript was used as exemplar or not is 
recorded in the `exemplar` column. 

Transcripts are grouped into genes based on the parameters 
`--p_exon_overlap`, and `--p_exons_overlap`. In order for two transcripts to 
be assigned to the same gene, a specified proportion of their exons must 
overlap. The proportional extent of overlap required between two exons is 
set using `--p_exon_overlap`. The proportion is calculated as the overlap 
divided by the shorter exon's length. The proportion of exons required to 
overlap to that degree is specified using `--p_exons_overlap`. That 
proportion is calculated as the number of overlapping exons divided by the 
number of exons in the transcript with fewer exons. For 
example, if `--p_exon_overlap` is set to `0.5`, then if one exon is 20 bases 
long and the other is 30 bases long, the overlap would need to be at least 
`0.5 * min(20, 30))`, or 10 bases long for the exons to have sufficient 
overlap to be counted. If we set `--p_exons_overlap` to `0.2`, then comparing 
a 10 exon transcript with a 15 exon transcript, at least `0.2 * min(10, 15)`, 
or 2 exons would need to sufficiently overlap for the transcripts to be 
assigned the same gene ID. It doesn't matter which exons match each other. If 
the third and sixth exons of the first transcript overlapped with the eighth 
and second exons, respectively, of the second transcript, the two transcripts 
would be assigned the same gene ID. For comparisons involving transcripts 
with 5 exons or fewer, at least one exon would need to overlap.

Transcripts in the output GTF file are assigned new potentially shared gene 
IDs and unique transcript IDs. By default, gene IDs take the form `LOC.M`, 
where `M` is a sequentially assigned integer (e.g. `LOC.1374`). Transcript IDs 
take the form `LOC.M.N`, where `LOC.M` corresponds to the associated gene ID, 
and `N` is an integer sequentially assigned to transcripts associated with 
that gene ID (e.g. `LOC.1374.2`). The prefix `LOC` can be changed using 
the `--gene_prefix` parameter. Mappings from original identifiers to new 
identifiers are found in `union.xrefs.tsv`. Transcripts that have unexpected 
structures (e.g. due to interchromosomal fusions or local genomic 
rearrangements) are identified based on the following conditions:

1) Exons on different chromosomes
2) Exons on different strands
3) Exon start coordinates out of order
4) Too much space between exons

The amount of space allowed between exons can be adjusted with the 
`--max_intron_length` parameter. Fusions are named with the prefix `fusion`, 
followed by a `.` delimited list of the genes to which the fusion maps to 
(e.g. `fusion.gene_id1.gene_id2`). Transcript rearrangement is flagged in the 
`rearranged` column of the `union.xrefs.tsv` file.



