import click
import pysam
import re
import os
import pandas as pd
from pathlib import Path
from typing import List
from dataclasses import dataclass
from .methylationPattern import methyl_patterns


@click.command()
@click.option('--inpath', type=click.Path(exists=True, readable=True),
              required=True, help='Directory with CpG and alignment files files')
@click.option('--thr', type=click.FLOAT, required=True, help='Threshold for '
                        'allele frequency bellow which an allele is ignored')
@click.option('--outpath', type=click.Path(writable=True), required=True,
              help='Path to place the output')
@click.option('--ampltable', type=click.Path(exists=True, readable=True),
              required=True, help="Tab separated file with amplicon locations")
def main(inpath, thr, outpath, ampltable):
    in_path = Path(inpath)
    alignment_files = list(in_path.glob("*bam"))

    # Loop over samples, put the per sample output in a dict according to
    # the amplicon
    ampl_to_df = {}
    for file in alignment_files:
        sample_id = sample_name(str(file))
        cpg_file = str(list(in_path.glob("CpG_OB_" + sample_id + "_bismark_*"))[0])
       # cpg_file = inpath + "/CpG_OB_" + sample_id + "_bismark_bt2.sorted.txt.gz"
        df = per_sample(file, thr, outpath, cpg_file, ampltable, sample_id)
        for ampl, d in df.items():
            if ampl not in ampl_to_df:
                ampl_to_df[ampl] = [d]
            else:
                ampl_to_df[ampl].append(d)

    # One table per amplicon
    for ampl, d in ampl_to_df.items():
        pd.concat(d).to_csv("{0}/{1}.tsv".format(outpath, ampl), sep ="\t",
                            header=True)
        pd.concat(d).to_excel("{0}/{1}.xls".format(outpath, ampl))
    # Create an empty file to signal the end of script execution for snakemake
    Path(outpath + '/methylator.txt').touch()


def per_sample(samfile, thr, outpath, cpgfile, ampltable, sample_id):
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

    # Create output directory
    if not os.path.exists(outpath):
        os.makedirs(outpath)

    # Load methylation call data. Forward and reverse strand are in two
    # separate files (OB and OT).
    # Combine them in one df. If file does not exist, create empty DF
    if os.path.isfile(cpgfile):
        methylationOB = pd.read_csv(cpgfile, sep="\t", skiprows=1,
                                    header=None, names=["Read", "MethylStatus", "Chr", "Pos", "Zz"])
    else:
        methylationOB = pd.DataFrame(columns=["Read", "MethylStatus", "Chr", "Pos", "Zz"])
    cpgfileOT = cpgfile.replace("CpG_OB_", "CpG_OT_")
    if os.path.isfile(cpgfileOT):
        methylationOT = pd.read_csv(cpgfileOT, sep="\t", skiprows=1,
                                    header=None, names=["Read", "MethylStatus", "Chr", "Pos", "Zz"])
        methylation = pd.concat([methylationOB, methylationOT])
    else:
        methylation = methylationOB

    # Read alignment file, create a dict of allele to records list
    samFile = pysam.AlignmentFile(samfile, 'rb')

    ampl_list = read_amplicon(ampltable)

    amplicon_to_df = {}

    # Loop over the list of amplicons
    for amplicon in ampl_list:
        start = amplicon.start
        end = amplicon.end
        methylation_thr = amplicon.methyl_thr
        amplicon_name = amplicon.name
        chrom = amplicon.chrom
        number_CGs = amplicon.nr_cg
        snp_coords = amplicon.snps_coord.split(";")
        # Extract the methylation table that concerns this amplicon region
        amplicon_meth = methylation.loc[(methylation["Chr"] == chrom) & (
                methylation["Pos"] >= start) & (methylation["Pos"] <= end)]

        index = []
        list_series = []

        # Check if amplicon has SNP
        if '-' not in amplicon.snps_coord:
            # Each amplicon region might have more than one SNPs
            for snp_coord in snp_coords: # TODO: what if there are no SNPs
                snp_coord = int(snp_coord)-1  # pysam uses zero-based indexing
                # Create a dict of allele to records list from the alignment file
                allele_to_read_record = base_to_reads(samFile, chrom, snp_coord)
                # returns a flattened list
                all_records = [record for list_records in
                               allele_to_read_record.values() for record in list_records]
                all_record_counts = len(all_records)

                # Remove alleles with read count below threshold
                records_to_keep = {}
                all_records_to_keep = []
                for allele, records in allele_to_read_record.items():
                    if len(records) > float(thr) * all_record_counts:
                        records_to_keep[allele] = records
                        all_records_to_keep.extend(records)

                # Add all alleles with all_records_to_keep to dict
                records_to_keep["Total"] = all_records_to_keep
                counts = phase_reads(records_to_keep, amplicon_meth,
                                     outpath, methylation_thr, number_CGs, sample_id, chrom, snp_coord)
                for allele, series in counts.items():
                    index.append((sample_id, amplicon_name, "{0}:{1}".format(chrom, snp_coord + 1), allele))
                    list_series.append(series)
            index = pd.MultiIndex.from_tuples(index, names=["Sample", "Amplicon", "SNP_coord", "Allele"])
            df = pd.DataFrame(list_series, index = index)

        else:
            series = methyl_patterns(amplicon_meth, outpath,
                                     methylation_thr, number_CGs, sample_id, "-", chrom, "-")
            index.append((sample_id, amplicon_name, "-", "-"))
            list_series.append(series)
            index = pd.MultiIndex.from_tuples(index, names=["Sample", "Amplicon", "SNP_coord", "Allele"])
            df = pd.DataFrame(list_series, index = index)

        if amplicon_name in amplicon_to_df:
            raise Exception("Amplicon name: {} is present twice. Check your ampltable and make "
                            "sure amplicon names are unique".format(amplicon_name))
        amplicon_to_df[amplicon_name] = df
    return(amplicon_to_df)


def base_to_reads(sam_file, chr, pos):
    """
    For all possible alleles in a genomic position, return a dictionary {
    allele => [Read IDs with this allele]}
    :param sam_file:
    :param chr:
    :param pos:
    :return: dictionary {allele => [Read IDs with this allele]}
    """
    pileups = sam_file.pileup(chr, pos, max_depth=30000)

    allele_to_read_record = {}
    for pileup_col in pileups:
        for pileup_read in pileup_col.pileups:
            if not pileup_read.is_del and not pileup_read.is_refskip and pileup_col.pos == pos:
                aln = pileup_read.alignment
                base = aln.query_sequence[pileup_read.query_position]
                if base not in allele_to_read_record:
                    allele_to_read_record[base] = [aln.query_name]
                else:
                    allele_to_read_record[base].append(aln.query_name)
    return(allele_to_read_record)


# A class to hold info about an amplicon
@dataclass
class Amplicon:
    name: str
    chrom: str
    start: int
    end: int
    strand: str
    nr_cg: int
    methyl_thr: int
    snps_coord: str


def read_amplicon(ampltable) -> List[Amplicon]:
    """
    Read amplicon table and get list of amplicon objects
    :param ampltable:
    :param methylation:
    :return:
    """

    # Load the amplicon table, loop over rows
    amplicons = pd.read_csv(ampltable, sep="\t")
    amplList = []
    for index, row in amplicons.iterrows():
        name = row["Name"]
        chrom = row["Chr"]
        start = row["start"]
        end = row["end"]
        strand = row["strand"]
        upper_mCG_thr = row["upper_mCG_thr"]
        low_mCG_thr = row["low_mCG_thr"]
        snp_coord = str(row["snps_coord"])
        # row: List[things]
        # amplication = Amplicon(*row)  <- Amplicon(row[0], row[1] .... )
        # row: Dict[str, things]
        # amplicon = Amplicon(**row) <- Amplicon(key=row[key], key2=row[key2]...)
        amplicon = Amplicon(name, chrom, start, end, strand, upper_mCG_thr, low_mCG_thr, snp_coord)
        amplList.append(amplicon)
    return amplList


def phase_reads(records_to_keep, methylation, outpath, low_mCG_thr, upper_mCG_thr, sample_id, chrom, snp_coord):
    """
    Provided a heterozygous SNP is present, phase reads according to SNP.
    Apply methylPatterns on split dataset
    :param outpath:
    :param cpgfile:
    :return: {allele => countsDF}
    """

    allele_to_counts = {}
    ## Loop over alleles, phase reads
    for allele, records in records_to_keep.items():
        methylation_phased = methylation[methylation["Read"].isin(records)]
        counts_per_class = methyl_patterns(methylation_phased, outpath,
                                           low_mCG_thr, upper_mCG_thr, sample_id, allele, chrom, snp_coord)
        allele_to_counts[allele] = counts_per_class

    return allele_to_counts


def sample_name(file):
    """
    Get sample name out of the bismark file name.
    Expects full path of a CpG file created by bismark
    """
    # str_search = re.search('.+/CpG_OB_(.+)_bismark.+', file)
    str_search = re.search('.+/(.+)_bismark_(bt2|hisat2)\.sorted\.bam', file)
    sample_name = str_search.group(1)
    return sample_name


if __name__ == '__main__':
    main()
