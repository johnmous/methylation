## Differential methylation analysis tool

A tool for the analysis of amplicon methylation data. Works with the output 
from [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/), a program to map bisulfite treated sequencing reads to a genome of interest and perform methylation calls.

Aims to:
  * Classify the reads accodring to their methylation status (number of CpGs that 
  are methylated) into: mostly methylated, mostly unmethylated and partially 
  methylated
  * Phase the reads according to a heterozygous SNP in the amplicon region
  * Count the number of reads in each methylation status class and report it per 
  group of phased reads
  
### Usage
Options:

| Option | Type | Description |
|:-:|:-:|:-:|
| --inpath | PATH | Directory with CpG and alignment files files  [required] |
| --thr | FLOAT | Threshold for allele frequency bellow which an allele is ignored [required]|
| --outpath | PATH | Path to place the output  [required] |
| --ampltable | PATH | Tab separated file with amplicon locations  [required] |

* **--inpath** a directory is provided where the output the alignment files 
and the methylation calls is placed. Script will first try to match files with glob
(*bam), extract the sample name and then look for the methylation call files
 (CpG*). All file names are expected to have the format as output by bismark 
 aligner and bismark_methylation_extractor
 
* **--ampltable** a tab separated file is provided, with the fields as shown
in the following example: 

| Name | Chr | start | end	| strand | size_bp | nr_CG | methyl_thr | snps_coord |
|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|:-:|
| H19 | chr11 |	1999754| 2000063 | neg | 310| 23 | 3 | 1999845;1999934 |

Most of the fields are self explanatory. The rest are as:
* **nr_CG** is an integer, the number of CpGs in the amplicon. 
* **methyl_thr** is an integer used to classify the reads as follows: If the 
number of 
methylated CpGs >= nr_CG - methyl_thr, then the read is methylated. If the 
number of methylated CpGs <=  methyl_thr the read is unmethylated. The rest 
are classified as partially methylated.
* **snps_coord** is an integer used to provide the location(s) of the 
heterozygous SNP(s) on which the reads are phased. If more than one, they 
should be separated with **;** 
