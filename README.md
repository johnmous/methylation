## Methylation analysis

* A collection of scripts for the analysis of methylation data
* Works with the output from [Bismark](https://www.bioinformatics.babraham.ac.uk/projects/bismark/), a program to map bisulfite treated sequencing reads to a genome of interest and perform methylation calls
* Meant for the analysis of amplicon data
* Aims to:
  * Classify the reads accodring to their methylation status (number of CpGs that 
  are methylated) into: mostly methylated, mostly unmethylated and partially 
  methylated
  * Phase the reads according to a heterozygous SNP in the amplicon region
  * Count the number of reads in each methylation status class and report it per 
  group of phased reads
