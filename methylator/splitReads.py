import click
import pysam
import re
import os
import pandas as pd
from pathlib import Path
from typing import List
from dataclasses import dataclass
import matplotlib
matplotlib.use('PDF') 
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from methylationPattern import methyl_patterns


@dataclass
class Per_sample_data:
    """
    Save the per sample data:
    A list with methylation data frames
    A list with Plot_data objects
    """
    methylation_list: list
    plot_list: list

@dataclass
class Plot_data:
    """
    Save data for plots
    """
    count_methyl_CpGs: pd.DataFrame
    sample_id: str
    amplicon: str
    chrom: str
    snp_coord: int
    low_mCG_thr: int
    upper_mCG_thr: int

@dataclass
class Methylation_data:
    amplicon: str
    data_frame: pd.DataFrame

# A class to hold info about an amplicon
@dataclass
class Amplicon:
    name: str
    chrom: str
    start: int
    end: int
    strand: str
    upper_mCG_thr: int
    low_mCG_thr: int
    snps_coord: str

@click.command()
@click.option('--inpath', type=click.Path(exists=True, readable=True),
              required=True, help='Directory with CpG and alignment files files')
@click.option('--thr', type=click.FLOAT, required=True, help='Threshold for '
                        'allele frequency bellow which an allele is ignored')
@click.option('--outpath', type=click.Path(writable=True), required=True,
              help='Path to place the output')
@click.option('--ampltable', type=click.Path(exists=True, readable=True),
              required=True, help="Tab separated file with amplicon locations")
@click.option('--plotgrid', type = str, default="3;2", help = 'Number of '
              'plots to draw per row and column, separated with ";"')
def main(inpath, thr, outpath, ampltable, plotgrid):
    in_path = Path(inpath)
    alignment_files = list(in_path.glob("*bam"))

    # Loop over samples
    ampl_to_df = {}
    ampl_snp_to_plot_data = {}
    for file in alignment_files:
        sample_id = sample_name(str(file))
        per_sample_data = per_sample(file, thr, in_path, outpath, ampltable,
                              sample_id)

        # Put the methylation tables in a dictionary, with amplicon name as
        # key
        methylation_data = per_sample_data.methylation_list
        for plot_data in methylation_data:
            if plot_data.amplicon not in ampl_to_df:
                ampl_to_df[plot_data.amplicon] = [plot_data.data_frame]
            else:
                ampl_to_df[plot_data.amplicon].append(plot_data.data_frame)

        # Put the plot data in a dictionary, amplicon name as key
        plot_data_list = per_sample_data.plot_list
        for plot_data in plot_data_list:
            snp_to_plot_data = {}
            ampl_snp = "{}.{}:{}".format(plot_data.amplicon, plot_data.chrom, \
                                                     str(plot_data.snp_coord))
            if ampl_snp not in ampl_snp_to_plot_data:
                ampl_snp_to_plot_data[ampl_snp] = [plot_data]
            else:
                ampl_snp_to_plot_data[ampl_snp].append(plot_data)

    # Write tables, per amplicon
    for ampl, d in ampl_to_df.items():
        pd.concat(d).to_csv("{0}/{1}.tsv".format(outpath, ampl), sep ="\t",
                            header=True)
        pd.concat(d).to_excel("{0}/{1}.xls".format(outpath, ampl))

    #
    plot_rows_per_page = int(plotgrid.split(";")[0])
    plot_columns_per_page = int(plotgrid.split(";")[1])
    plots_per_page = plot_rows_per_page*plot_columns_per_page

    # Make all plots per amplicon and snp and save in one PDF
    for ampl_snp, plot_data_list in ampl_snp_to_plot_data.items():
        # Make plots
        matplotlib.rcParams.update({'font.size': 5})
        pdf =  PdfPages("{}/plots/{}.pdf".format(outpath, ampl_snp))
        fig = plt.figure()
        plot_count = 0
        for plot_data in plot_data_list:
            count_methyl_CpcGs =  plot_data.count_methyl_CpGs
            low_mCG_thr = plot_data.low_mCG_thr
            upper_mCG_thr = plot_data.upper_mCG_thr
            amplicon = plot_data.amplicon
            chrom = plot_data.chrom
            snp_coord = plot_data.snp_coord
            sample_id = plot_data.sample_id
            if len(count_methyl_CpcGs.index) > 0:
                plot_count += 1
                ax = fig.add_subplot(plot_rows_per_page,
                                     plot_columns_per_page, plot_count)
                ax.plot("methStatesCount", "Total", data=count_methyl_CpcGs,
                        label="Total", color="black", linestyle=":")
                alleles = count_methyl_CpcGs.drop(
                    labels=["methStatesCount", "Total"], axis=1)

                # Each allele has a fixed color
                colors = {"C": "blue",
                          "T": "red",
                          "G": "yellow",
                          "A": "green"}
                for allele in list(alleles.columns):
                    ax.scatter(x=count_methyl_CpcGs["methStatesCount"],
                                y=count_methyl_CpcGs[allele], label=allele,
                                color=colors[allele], s = 10)
                ax.legend()
                ax.grid(True)
                ax.axvline(low_mCG_thr, color="black", linestyle="--")
                ax.axvline(upper_mCG_thr, color="black", linestyle="--")
                plt.title("{}_{}_{}:{}".format(amplicon, sample_id, chrom,
                                               snp_coord))
                plt.xlim(left=-1, right=25)
                plt.xlabel("Number of methylated CpG sites")
                plt.ylabel("Number of reads")

            # If max number of plots per page reached, save figure (page)
            # and create new figure object
            if plot_count == plots_per_page:
                fig.tight_layout()
                pdf.savefig(fig)
                fig = plt.figure()
                plot_count = 0
        fig.tight_layout()
        pdf.savefig(fig)
        pdf.close()
    # Create an empty file to signal the end of script execution for snakemake
    Path(outpath + '/methylator.txt').touch()

def per_sample(samfile, thr, in_path, outpath, ampltable, sample_id):
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
    if not os.path.exists(outpath + "/plots"):
        os.makedirs(outpath+ "/plots")

    # Load methylation call data. Forward and reverse strand are in two
    # separate files (OB and OT).
    # Combine them in one df. If file does not exist, create empty DF
    cpg_OB_file = list(in_path.glob("CpG_OB*" + sample_id + "_bismark_*"))
    cpg_OT_file = list(in_path.glob("CpG_OT*" + sample_id + "_bismark_*"))
    # List should contain one file exactly
    if len(cpg_OB_file) == 1:
        cpg_string = str(cpg_OB_file[0])
        methylationOB = pd.read_csv(cpg_string, sep="\t", skiprows=1,
                                    header=None, names=["Read", "MethylStatus", "Chr", "Pos", "Zz"])
    else:
        methylationOB = pd.DataFrame(columns=["Read", "MethylStatus", "Chr", "Pos", "Zz"])
    if len(cpg_OT_file) == 1:
        cpg_string = str(cpg_OT_file[0])
        methylationOT = pd.read_csv(cpg_string, sep="\t", skiprows=1,
                                    header=None, names=["Read", "MethylStatus", "Chr", "Pos", "Zz"])
        methylation = pd.concat([methylationOB, methylationOT])
    else:
        methylation = methylationOB

    # Read alignment file, create a dict of allele to records list
    samFile = pysam.AlignmentFile(samfile, 'rb')
    ampl_list = read_amplicon(ampltable)

    methyl_data_list = []
    plot_data_list = []
    # Loop over the list of amplicons
    for amplicon in ampl_list:
        start = amplicon.start
        end = amplicon.end
        low_mCG_thr = amplicon.low_mCG_thr
        amplicon_name = amplicon.name
        chrom = amplicon.chrom
        upper_mCG_thr = amplicon.upper_mCG_thr
        snp_coords = amplicon.snps_coord.split(";")

        # Extract the methylation table that concerns this amplicon region
        amplicon_meth = methylation.loc[(methylation["Chr"] == chrom) & (
                methylation["Pos"] >= start) & (methylation["Pos"] <= end)]

        index = []
        list_series = []
        # Check if amplicon has SNP
        if '-' not in amplicon.snps_coord:
            # Each amplicon region might have more than one SNPs
            for snp_coord in snp_coords:
                snp_coord = int(snp_coord)-1  # pysam uses zero-based indexing
                # Create a dict of allele (base) to a list of records from
                # the alignment file
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

                # Loop over alleles
                allele_to_counts = {}
                data_frames_plots = []
                # Loop over alleles, phase reads
                for allele, records in records_to_keep.items():
                    methylation_phased = methylation[methylation["Read"].isin(records)]
                    methyl_DFs = methyl_patterns(methylation_phased, outpath,
                                                 low_mCG_thr, upper_mCG_thr,
                                                 sample_id, allele, chrom,
                                                 snp_coord)
                    allele_to_counts[allele] = methyl_DFs.count_meth_class
                    data_frames_plots.append(methyl_DFs.count_methyl_CpGs)
                count_methyl_CpGs = pd.concat(data_frames_plots, axis=1,
                                              join="outer", keys="methStatesCount")
                count_methyl_CpGs.columns = count_methyl_CpGs.columns.droplevel()
                count_methyl_CpGs.reset_index(inplace=True)
                count_methyl_CpGs.fillna(value=0, inplace=True)
                plot_data = Plot_data(count_methyl_CpGs, sample_id, amplicon_name,
                                      chrom, snp_coord, low_mCG_thr,
                                      upper_mCG_thr)
                plot_data_list.append(plot_data)

                for allele, series in allele_to_counts.items():
                    index.append((sample_id, amplicon_name, "{0}:{1}".format(chrom, snp_coord + 1), allele))
                    list_series.append(series)
            index = pd.MultiIndex.from_tuples(index, names=["Sample", "Amplicon", "SNP_coord", "Allele"])
            df = pd.DataFrame(list_series, index = index)

        else:
            series = methyl_patterns(amplicon_meth, outpath,
                                     low_mCG_thr, upper_mCG_thr, sample_id, "-", chrom, "-")
            index.append((sample_id, amplicon_name, "-", "-"))
            list_series.append(series)
            index = pd.MultiIndex.from_tuples(index, names=["Sample", "Amplicon", "SNP_coord", "Allele"])
            df = pd.DataFrame(list_series, index = index)
            # TODO: plot_data

        # if amplicon_name in amplicon_to_df:
        #    raise Exception("Amplicon name: {} is present twice. Check your ampltable and make "
        #                    "sure amplicon names are unique".format(
        #                    amplicon_name))
        # amplicon_to_df[amplicon_name] = df
        methylation_data = Methylation_data(amplicon_name, df)
        methyl_data_list.append(methylation_data)
    per_sample_data = Per_sample_data(methyl_data_list, plot_data_list)
    return(per_sample_data)


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


def sample_name(file):
    """
    Get sample name out of the bismark file name.
    Expects full path of a CpG file created by bismark
    """
    str_search = re.search('.+/(.+)_bismark_(bt2|hisat2)\.sorted\.bam', file)
    sample_name = str_search.group(1)
    return sample_name


if __name__ == '__main__':
    main()
