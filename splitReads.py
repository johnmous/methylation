## Title: Group reads based on heterozygous SNPs locations. Each group of reads will be placed in a separate BAM file
## Author: I. Moustakas, i.moustakas@lumc.nl

import click
import pysam
import re
import os
import pandas as pd
from methylationPattern import methylPatterns

@click.command()
@click.option('--samfile', required=True, help='Alignment file')
@click.option('--thr', required=True, help='Threshold for allele frequency bellow which an allele is ignored')
@click.option('--outpath', required=True, help='Path to place the output')
@click.option('--cpgfile', required=True, help="CpG file from bismark (CpG_OB_*)")
@click.option('--ampltable', required=True, help="Tab separated file with amplicon locations") ##TODO: specify table format
def perSample(samfile, thr, outpath, cpgfile, ampltable):
    """
    Accepts an alignment file, a methylation call file and list of amplicons
    For each amplicon, calculate a table with methylation pattern counts
    :param samfiles:
    :param location:
    :param thr:
    :param outpath:
    :param cpgfile:
    :param ampltable:
    :return: {amplicon => methylation counts DF}
    """

    ## Create output directory
    if not os.path.exists(outpath):
        os.makedirs(outpath)

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

    #amplInfo = getAmplicon(ampltable, methylation)
    amplList = readAmplicon(ampltable)
    sampleNm = sampleName(cpgfile)

    ampliconToDF = {}

    for amplicon in amplList:
        start = amplicon.start
        end = amplicon.end
        methylationThr = amplicon.methylThr
        ampliconName = amplicon.name
        chrom = amplicon.chrom
        snpCoords = amplicon.snp_coord.split(";")

        ## Extract the methylation table that concerns this amplicon region
        ampliconMeth = methylation.loc[(methylation["Chr"] == chrom) & (methylation["Pos"] >= start) & (methylation["Pos"] <= end)]

        index = []
        listSeries = []
        ## Each amplicon region might have more than one SNPs
        for snpCoord in snpCoords:
            snpCoord = int(snpCoord)-1
            ## Create a dict of allele to records list from the alignment file
            alleleToReadRecord = baseToReads(samFile, chrom, snpCoord)
            allRecords = [record for listRecords in alleleToReadRecord.values() for record in listRecords]
            allRecordCounts = len(allRecords)

            ## Remove alleles with read count below threshold
            recordsToKeep = {}
            for allele, records in alleleToReadRecord.items():
                if len(records) > float(thr) * allRecordCounts:
                    recordsToKeep[allele] = records

            ## Add All Alleles with allRecords to dict
            recordsToKeep["All"] = allRecords

            counts = phaseReads(recordsToKeep, ampliconMeth, outpath, methylationThr)
            for allele, series in counts.items():
                index.append((sampleNm, ampliconName, "{0}:{1}".format(chrom, snpCoord), allele))
                listSeries.append(series)
        index = pd.MultiIndex.from_tuples(index, names=["Sample", "Amplicon", "SNP_coord", "Allele"])
        df = pd.DataFrame(listSeries, index = index)
        ampliconToDF[ampliconName] = df

    print(df)
    return(ampliconToDF)


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


## A class to hold info about an amplicon
class Amplicon(object):
    def __init__(self, name, chrom, start, end, strand, methylThr, snp_coord):
        self.name = name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self.methylThr = methylThr
        self.snp_coord = snp_coord


def readAmplicon(ampltable):
    """
    Read amplicon table and get list of amplicon objects
    :param ampltable:
    :param methylation:
    :return:
    """
    ## Load the amplicon table
    amplicons = pd.read_csv(ampltable, sep="\t")

    amplList = []
    for index, row in amplicons.iterrows():
        name = row["Name"]
        chrom = row["Chr"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]
        methylThr = row["methylThr"]
        snp_coord = row["snp_coord"]
        amplicon = Amplicon(name, chrom, start, end, strand, methylThr, snp_coord)
        amplList.append(amplicon)


    return amplList


def SplitBismMethExtr(methylExtr, readIDs):
    """
    Split Bismark methylation extractor output file.
    Phase the reads according to the SNP(s) in the amplicon and
    save in separate file
    """


def phaseReads(recordsToKeep, methylation, outpath, methylThr):
    """
    Provided a heterozygous SNP is present, phase reads according to SNP.
    Apply methylPatterns on split dataset
    :param outpath:
    :param cpgfile:
    :return: {allele => countsDF}
    """

    alleleToCounts = {}
    ## Loop over alleles, phase reads
    for allele, records in recordsToKeep.items():
        methylationPhased = methylation[methylation["Read"].isin(records)]
        countsPerClass = methylPatterns(methylationPhased, outpath, methylThr)
        alleleToCounts[allele] = countsPerClass

    return alleleToCounts

def sampleName(file):
    """
    Get sample name out of the bismark file name.
    Expects full path of a CpG file created by bismark
    """
    str_search = re.search('.+/CpG_OB_(.+)_bismark.+', file)
    sampleName = str_search.group(1)
    return sampleName

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

