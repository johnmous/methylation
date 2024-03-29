# Copyright (c) 2022 Leiden University Medical Center
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

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
from .methylationPattern import methyl_patterns


@dataclass
class Per_sample_data:
    """
    Save the per sample data:
    A list with methylation data frames
    A list with Plot_data objects
    """
    methylation_list: list
    plot_list: list
    positional_meth_pct_table_per_sample: pd.DataFrame


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
    ampl_prev_thr: float


@click.command()
@click.option('--inpath', type=click.Path(exists=True, readable=True),
              required=True, help='Directory with CpG and alignment files files')
@click.option('--thr', type=click.FLOAT, required=True,
              help='Threshold for allele frequency bellow which an '
                   'allele is ignored')
@click.option('--outpath', type=click.Path(writable=True), required=True,
              help='Path to place the output')
@click.option('--ampltable', type=click.Path(exists=True, readable=True),
              required=True, help="Tab separated file with amplicon "
                                  "locations on genome")
@click.option('--plotgrid', type=str, default="3;2",
              help='Number of plots to draw per row and column, '
                   'separated with ";"')
def main(inpath, thr, outpath, ampltable, plotgrid):
    in_path = Path(inpath)
    alignment_files = list(in_path.glob("*bam"))

    # Initiate dictionaries and DF to store info and loop over samples
    ampl_to_df = {}
    ampl_snp_to_plot_data = {}
    positional_meth_pct_table = pd.DataFrame()
    for file in alignment_files:
        sample_id = sample_name(str(file))
        print(f"Working with sample: {sample_id}")
        per_sample_data = per_sample(file, thr, in_path, outpath, ampltable,
                                     sample_id)

        # Positional methylation percentages
        positional_meth_pct_table_per_sample = per_sample_data.\
            positional_meth_pct_table_per_sample
        positional_meth_pct_table_per_sample["sample"] = sample_id
        positional_meth_pct_table = pd.concat([positional_meth_pct_table,
                                               positional_meth_pct_table_per_sample])

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
            ampl_snp = f"{plot_data.amplicon}.{plot_data.chrom}:{str(plot_data.snp_coord)}"
            if ampl_snp not in ampl_snp_to_plot_data:
                ampl_snp_to_plot_data[ampl_snp] = [plot_data]
            else:
                ampl_snp_to_plot_data[ampl_snp].append(plot_data)

    positional_meth_pct_table.columns.name = None
    positional_meth_pct_table.to_csv('positional_meth_pct_table.csv')

    # Write tables, per amplicon
    for ampl, d in ampl_to_df.items():
        pd.concat(d).to_csv(f"{outpath}/{ampl}.tsv", sep="\t",
                            header=True)
        pd.concat(d).to_excel(f"{outpath}/{ampl}.xls")

    #
    plot_rows_per_page = int(plotgrid.split(";")[0])
    plot_columns_per_page = int(plotgrid.split(";")[1])
    plots_per_page = plot_rows_per_page * plot_columns_per_page

    # Make all plots per amplicon and snp and save in one PDF
    for ampl_snp, plot_data_list in ampl_snp_to_plot_data.items():
        matplotlib.rcParams.update({'font.size': 5})
        pdf = PdfPages(f"{outpath}/plots/{ampl_snp}.pdf")
        fig = plt.figure()
        plot_count = 0
        for plot_data in plot_data_list:
            count_methyl_CpcGs = plot_data.count_methyl_CpGs
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
                          "A": "green",
                          "-": "black"}
                for allele in list(alleles.columns):
                    ax.scatter(x=count_methyl_CpcGs["methStatesCount"],
                               y=count_methyl_CpcGs[allele], label=allele,
                               color=colors[allele], s=10)
                ax.legend()
                ax.grid(True)
                ax.axvline(low_mCG_thr, color="black", linestyle="--")
                ax.axvline(upper_mCG_thr, color="black", linestyle="--")
                plt.title(f"{amplicon}_{sample_id}_{chrom}:{snp_coord}")
                plt.xlim(left=-1, right=25)
                plt.xlabel("Number of per molecule methylated CpG sites")
                plt.ylabel("Number of molecules")

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

    # Plot the per position methylation percentages
    pdf = PdfPages(f"{outpath}/plots/{amplicon}_methylation_pct.pdf")
    positional_meth_pct_table = positional_meth_pct_table[
        positional_meth_pct_table['allele'] == "Total"]
    positional_meth_pct_table.sort_values('sample', inplace=True)
    samples = positional_meth_pct_table["sample"]
    positional_meth_pct_table = positional_meth_pct_table.\
        drop(["allele", "amplicon", "sample"], axis=1)
    positional_meth_pct_table = positional_meth_pct_table * 100
    positions = positional_meth_pct_table.columns

    fig, ax = plt.subplots()
    im = ax.imshow(positional_meth_pct_table)

    # Show all ticks and label them with the respective list entries
    ax.set_yticks(range(0, len(samples)))
    ax.set_yticklabels(samples)
    ax.set_xticks(range(0, len(positions)))
    ax.set_xticklabels(positions)
    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("% methylated CpGs", rotation=-90, va="bottom")

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=90, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(positional_meth_pct_table.shape[0]):
        for j in range(positional_meth_pct_table.shape[1]):
            text = ax.text(j, i, int(positional_meth_pct_table.iloc[i, j]),
                           ha="center", va="center", color="r", fontsize=2.5)
    ax.set_title("Per position % of methylated CpGs")
    fig.tight_layout()
    pdf.savefig(fig)
    pdf.close()

    # Create an empty file to signal the end of script execution for snakemake
    Path(outpath + '/methylator.txt').touch()


def per_sample(sam_file, thr, in_path, out_path, ampl_table, sample_id):
    """
    Accepts an alignment file, a methylation call file and list of amplicons
    For each amplicon, calculate a table with methylation pattern counts
    :param sam_file: Alignment file
    :param thr: Threshold for allele frequency bellow which an allele is
    ignored
    :param in_path: Input path where alignment and CpG count files are
    located
    :param out_path: Path to store output
    :param ampl_table: Path of table with amplicon info
    :param sample_id: Sample name
    :return: {amplicon => methylation counts DF}
    """

    # Create output directory
    if not os.path.exists(out_path):
        os.makedirs(out_path)
    if not os.path.exists(out_path + "/plots"):
        os.makedirs(out_path + "/plots")

    # Load methylation call data. Forward and reverse strand are in two
    # separate files (OB and OT).
    # Combine them in one df. If file does not exist, create empty DF
    cpg_OB_file = list(in_path.glob("CpG_OB*" + sample_id + "_bismark_*"))
    cpg_OT_file = list(in_path.glob("CpG_OT*" + sample_id + "_bismark_*"))
    # List should contain one file exactly
    if len(cpg_OB_file) == 1:
        cpg_string = str(cpg_OB_file[0])
        methylationOB = pd.read_csv(cpg_string, sep="\t", skiprows=1,
                                    header=None, names=["Read",
                                                        "MethylStatus", "Chr",
                                                        "Pos", "Zz"])
    else:
        methylationOB = pd.DataFrame(columns=["Read", "MethylStatus", "Chr",
                                              "Pos", "Zz"])
    if len(cpg_OT_file) == 1:
        cpg_string = str(cpg_OT_file[0])
        methylationOT = pd.read_csv(cpg_string, sep="\t", skiprows=1,
                                    header=None, names=["Read", "MethylStatus",
                                                        "Chr", "Pos", "Zz"])
        methylation = pd.concat([methylationOB, methylationOT])
    else:
        methylation = methylationOB

    # Read alignment file, create a dict of allele to records list
    samFile = pysam.AlignmentFile(sam_file, 'rb')
    ampl_list = read_amplicon(ampl_table)

    methyl_data_list = []
    plot_data_list = []
    positional_meth_pct_table_per_sample = pd.DataFrame()

    # Loop over the list of amplicons
    for amplicon in ampl_list:
        start = amplicon.start
        end = amplicon.end
        low_mCG_thr = amplicon.low_mCG_thr
        amplicon_name = amplicon.name
        chrom = amplicon.chrom
        upper_mCG_thr = amplicon.upper_mCG_thr
        snp_coords = amplicon.snps_coord.split(";")
        ampl_prev_thr = amplicon.ampl_prev_thr
        print(f"Starting with amplicon: {amplicon_name}")

        # Extract the methylation table that concerns this amplicon region
        amplicon_meth = methylation.loc[(methylation["Chr"] == chrom) & (
                methylation["Pos"] >= start) & (methylation["Pos"] <= end)]

        index = []
        list_series = []
        # Check if amplicon has SNP
        if '-' not in amplicon.snps_coord:
            # Each amplicon region might have more than one SNPs
            for snp_coord in snp_coords:
                snp_coord = int(snp_coord)
                # Create a dict of allele (base) to a list of records from
                # the alignment file
                allele_to_read_record = base_to_reads(samFile, chrom,
                                                      snp_coord - 1)  # pysam
                # uses zero-based indexing
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

                allele_to_counts = {}
                data_frames_plots = []
                # Loop over alleles, phase reads
                for allele, records in records_to_keep.items():
                    methylation_phased = amplicon_meth[amplicon_meth["Read"].isin(records)]
                    methyl_DFs = methyl_patterns(methylation_phased, out_path,
                                                 low_mCG_thr, upper_mCG_thr,
                                                 sample_id, allele, chrom,
                                                 snp_coord, ampl_prev_thr)
                    data_frames_plots.append(methyl_DFs.count_methyl_CpGs)
                    allele_to_counts[allele] = methyl_DFs.count_meth_class
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
                    index.append((sample_id, amplicon_name, f"{chrom}:{snp_coord}", allele))
                    list_series.append(series)
            index = pd.MultiIndex.from_tuples(index, names=["Sample",
                                                            "Amplicon",
                                                            "SNP_coord",
                                                            "Allele"])
            df = pd.DataFrame(list_series, index=index)

        # No SNP available
        else:
            # To improve performance, randomly select a number of records
            # (reads) to keep with the help of pysam
            position = int(abs((end - start) / 2))  # Middle of the amplicon
            print("start reading pileups")
            pileups = samFile.pileup(chrom, position, max_depth=3000)
            read_ids = []
            for pileup_col in pileups:
                for pileup_read in pileup_col.pileups:
                    if not pileup_read.is_del and not pileup_read.is_refskip:
                        name = pileup_read.alignment.query_name
                        if name not in read_ids:
                            read_ids.append(name)

            amplicon_meth = amplicon_meth[amplicon_meth["Read"].isin(read_ids)]

            # Plot data
            allele_to_counts = {}
            data_frames_plots = []
            snp_coord = "-"
            # Total column present in the DF is necessary for plotting
            for allele in ["-", "Total"]:
                methyl_DFs = methyl_patterns(amplicon_meth, out_path,
                                             low_mCG_thr, upper_mCG_thr,
                                             sample_id, allele, chrom,
                                             snp_coord, ampl_prev_thr)
                data_frames_plots.append(methyl_DFs.count_methyl_CpGs)
                allele_to_counts[allele] = methyl_DFs.count_meth_class

                positional_meth_pct = methyl_DFs.positional_methyl_pct
                positional_meth_pct['amplicon'] = amplicon_name
                positional_meth_pct['allele'] = allele
                positional_meth_pct_table_per_sample = pd.\
                    concat([positional_meth_pct_table_per_sample, positional_meth_pct])
            count_methyl_CpGs = pd.concat(data_frames_plots, axis=1,
                                          join="outer", keys="methStatesCount")
            count_methyl_CpGs.columns = count_methyl_CpGs.columns.droplevel()
            count_methyl_CpGs.reset_index(inplace=True)
            count_methyl_CpGs.fillna(value=0, inplace=True)
            plot_data = Plot_data(count_methyl_CpGs, sample_id, amplicon_name,
                                  chrom, snp_coord, low_mCG_thr,
                                  upper_mCG_thr)
            plot_data_list = [plot_data]
            for allele, series in allele_to_counts.items():
                index.append((sample_id, amplicon_name, f"{chrom}:{snp_coord}", allele))
                list_series.append(series)
            index = pd.MultiIndex.from_tuples(index, names=["Sample",
                                                            "Amplicon",
                                                            "SNP_coord",
                                                            "Allele"])
            df = pd.DataFrame(list_series, index=index)
        methylation_data = Methylation_data(amplicon_name, df)
        methyl_data_list.append(methylation_data)
    per_sample_data = Per_sample_data(methyl_data_list, plot_data_list,
                                      positional_meth_pct_table_per_sample)
    return per_sample_data


def base_to_reads(sam_file, chrom, pos):
    """
    For all possible alleles in a genomic position, return a dictionary {
    allele => [Read IDs with this allele]}
    :param sam_file: Alignment file
    :param chrom: Chromosome where SNP resides
    :param pos: Genomic position of SNP
    :return: dictionary {allele => [Read IDs with this allele]}
    """
    pileups = sam_file.pileup(chrom, pos, max_depth=3000)

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
    return (allele_to_read_record)


def read_amplicon(ampl_table) -> List[Amplicon]:
    """
    Read amplicon table and get list of amplicon objects
    :param ampl_table: Path of table with amplicon info
    :return: List[Amplicon]
    """

    # Load the amplicon table, loop over rows
    amplicons = pd.read_csv(ampl_table, sep="\t")
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
        ampl_prev_thr = row["ampl_prev_thr"]
        amplicon = Amplicon(name, chrom, start, end, strand, upper_mCG_thr,
                            low_mCG_thr, snp_coord, ampl_prev_thr)
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
