## Differential methylation analysis tool

A tool for the analysis of amplicon methylation data. Works with the output 
from [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/), a program used to map bisulfite treated sequencing reads to the genome of interest and perform methylation calls.

The tool delivers the following:
  * Classify the genomic sequences by their methylation status into: mostly methylated, mostly unmethylated and partially 
  methylated. Methylation status refers to the number of CpGs in a genomic sequence that are methylated
  * Phase the genomic sequences according to a heterozygous SNP in the amplicon region
  * Count the number of genomic sequences in each methylation status class and report it per 
  group of phased reads
  
### Usage
Options:

| Option | Type | Description |
|:-:|:-:|:-:|
| --inpath | PATH | Directory with CpG and alignment files files  [required] |
| --thr | FLOAT | Threshold for allele frequency bellow which an allele is ignored [required]|
| --outpath | PATH | Path to place the output  [required] |
| --ampltable | PATH | Tab separated file with amplicon attributes  [required] |
| --plotgrid | string | Number of plots to draw per row and column, separated with ";" | 

* **--inpath** a directory is provided where the output the alignment files 
and the methylation calls is placed. Script will first try to match files with glob
(*bam), extract the sample name and then look for the methylation call files
 (CpG*). All file names are expected to have the format as output by bismark 
 aligner and bismark_methylation_extractor
 
* **--ampltable** a tab separated file is provided, with the fields as shown
in the following example: 

| Name | Chr | start | end	| strand | size_bp | upper_mCG_thr | low_mCG_thr | snps_coord |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| H19 | chr11 |	1999754 | 2000063 | neg | 310 | 20 | 3 | 1999845;1999934 |

Most of the fields are self explanatory. The rest are as:
* **upper_mCG_CG** is an integer. If a read has this number, or more, CpGs 
methylated, it is classified as methylated. 
* **low_mCG_thr** is an integer. If a read has this number, or less, CpGs 
methylated, it is classified as unmethylated.
* **snps_coord** is an integer used to provide the location(s) of the 
heterozygous SNP(s) on which the reads are phased. If more than one, they 
should be separated with **;** 
