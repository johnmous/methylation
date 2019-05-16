#!/exports/sasc/ioannis/miniconda3/envs/snake/bin/python
## Title: Group reads based on heterozygous SNPs locations. Each group of reads will be placed in a separate BAM file
## Author: I. Moustakas, i.moustakas@lumc.nl

import click
import pysam
import re
import os
import subprocess
import pandas as pd
from methylationPattern import methylPatterns

@click.command()
@click.option('--samfile', required = True, help='Alignment file')
@click.option('--location', required = True, help='Allele location to split on, formatted as chrom:pos')
@click.option('--thr', required = True, help='Threshold for allele frequency bellow which this allele is ignored')
@click.option('--outpath', required = True, help='Path to place the output')
@click.option('--cpgfile', required = True, help = "CpG file from bismark (CpG_OB_*)")
@click.option('--ampltable', required = True, help = "Tab separated file with amplicon locations") ##TODO: specify table format
def perSample(samfile, location, thr, outpath, cpgfile, ampltable):
    ## Get command line arguments
    chrom, pos = location.split(":")
    pos = int(pos) -1

    ## Load the amplicon table
    amplicons = pd.read_csv(ampltable, sep="\t")

    ## Load methylation call data. Forward and reverse strand are in two separate files (OB and OT).
    ## Combine them in one df. If file does not exist, create empty DF
    if os.path.isfile(cpgfile):
        methylationOB = pd.read_csv(cpgfile, sep="\t", skiprows=1, header=None, names=["Read", "MethylStatus", "Chr", "Pos", "Zz"])
    else:
        methylationOB = pd.DataFrame(columns=["Read", "MethylStatus", "Chr", "Pos", "Zz"])
    cpgfileOT = cpgfile.replace("CpG_OB_", "CpG_OT_")
    if os.path.isfile(cpgfileOT):
        methylationOT = pd.read_csv(cpgfileOT, sep="\t", skiprows=1, header=None, names=["Read", "MethylStatus", "Chr", "Pos", "Zz"])
        methylation = pd.concat([methylationOB, methylationOT])
    else:
        methylation = methylationOB

    ## Read alignment file, create a dict of allele to records list
    samFile = pysam.AlignmentFile(samfile, 'rb')
    alleleToReadRecord = baseToReads(samFile, chrom, pos)
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    allRecordCounts = len([record for listRecords in alleleToReadRecord.values() for record in listRecords])

    ## Remove alleles with read count below threshold
    recordsToKeep = {}
    for allele, records in alleleToReadRecord.items():
        if len(records) > float(thr)*allRecordCounts:
            recordsToKeep[allele] = records

    amplToMeth = getAmplicon(ampltable, methylation)

    for name, methylation in amplToMeth.items():
        counts = phaseReads(recordsToKeep, methylation, outpath)
        print(counts)


def baseToReads(samFile, chr, pos):
    """
    For all possible alleles in a genomic position, return a dictionary {allele => [Read IDs with this allele]}
    :param samFile:
    :param chr:
    :param pos:
    :return: dictionary {allele => [Read IDs with this allele]}
    """
    pileups = samFile.pileup(chr, pos, max_depth=30000)

    baseToReadRecord = {}
    for pileupCol in pileups:
        for pileupRead in pileupCol.pileups:
            if not pileupRead.is_del and not pileupRead.is_refskip and pileupCol.pos == pos:
                aln = pileupRead.alignment
                base = aln.query_sequence[pileupRead.query_position]
                if base not in baseToReadRecord:
                    baseToReadRecord[base] = [aln.query_name]
                else:
                    baseToReadRecord[base].append(aln.query_name)
    return(baseToReadRecord)


def getAmplicon(ampltable, methylation):
    """
    Split the methylation table per amplicon.
    Return a dictionary
    :param ampltable:
    :param methylation:
    :return: {amplicon name => methylation table}
    """
    ## Load the amplicon table
    amplicons = pd.read_csv(ampltable, sep="\t")

    amplToMeth = {}
    ## Go through amplicons and extract from meth table
    for index, row in amplicons.iterrows():
        name = row["Name"]
        chrom = row["Chr"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]
        methylThr = row["methylThr"]

        ## Extract the methylation table that concerns this amplicon region
        ampliconMeth = methylation.loc[(methylation["Chr"] == chrom) & (methylation["Pos"] >= start) & (methylation["Pos"] <= end)]
        amplToMeth[name] = ampliconMeth

    return amplToMeth


def SplitBismMethExtr(methylExtr, readIDs):
    """
    Split Bismark methylation extractor output file.
    Phase the reads according to the SNP(s) in the amplicon and
    save in separate file
    """


def phaseReads(recordsToKeep, methylation, outpath):
    """
    Provided a heterozygous SNP is present, phase reads according to SNP.
    Apply methylPatterns on split dataset
    :param outpath:
    :param cpgfile:
    :return: {allele => countsDF}
    """

    ## Get the sample name out of the bam file name. Expects a Bismark output file
    ## str_search = re.search('.*/(.+)_bismark_bt2.sorted.bam', alignFile)
    ## sampleName = str_search.group(1)

    alleleToCounts = {}
    ## Loop over alleles, phase reads
    for allele, records in recordsToKeep.items():
        methylationPhased = methylation[methylation["Read"].isin(records)]
        countsPerClass = methylPatterns(methylationPhased, outpath)
        alleleToCounts[allele] = countsPerClass

    return alleleToCounts


if __name__ == '__main__':
    perSample()

####### $$$$$$$ #######
        # alleleSpecificFileID = "{0}_{1}_{2}_{3}".format(sampleName, chr, pos, allele)
        #
        # ## outputBam = "{0}/{1}_bismark_bt2.bam".format(outputPath, alleleSpecificFileID)
        # ## outputSortedBam = "{0}/bamFiles/{1}_bismark_bt2.sorted.bam".format(outputPath, alleleSpecificFileID)
        # if not os.path.exists(outputPath+"/bamFiles"):
        #
        #     os.makedirs(outputPath+"/bamFiles")
        # alleleSpecReadsBam = pysam.AlignmentFile(outputBam, "wb", template=samFile)
        # for rcd in recordsToKeep[allele]:
        #     alleleSpecReadsBam.write(rcd)
        # alleleSpecReadsBam.close()
        # pysam.sort("-o", outputSortedBam, outputBam)
        # pysam.index(outputSortedBam)
