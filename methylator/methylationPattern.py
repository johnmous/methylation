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

import pandas as pd
from pathlib import Path
from dataclasses import dataclass
import numpy as np

# Class to save the result of methyl_patterns
@dataclass
class MethylationDFs:
    count_meth_class: pd.Series
    count_methyl_CpGs: pd.DataFrame


def methyl_patterns(methyl_extr, outpath, methyl_thr, upper_mCG_thr, sample_id, allele, chrom, snp_coord, ampl_prev_thr):
    """
    Separate reads in three categories: Methylated , unmethylated, partially methylated according to the number
    of CpGs that are methyalted
    :param methyl_extr: pandas df with columns "Read", "MethylStatus", "Chr", "Pos", "Zz". This is the output of bismark_methylation_extractor
    :param outpath: path to save the table with the methylation patterns
    :param: number_CGs: Number of CGs in the amplicon
    :param methyl_thr: Integer to subtract from number_CGs to set threshold from mostly methylated
    :param: ampl_prev_thr: amplicon prevalence threshold, a float to remove CpG sites less than this fraction of total reads
    :return: A pandas.Series with read counts and percentages for the three methylation categories
    """

    # if records for amplicon
    if len(methyl_extr.index)>0:

        # Get all methylation positions in the amplicon
        meth_pos_counts = methyl_extr["Pos"].value_counts()
        read_pos_methyl = methyl_extr[["Read", "Pos", "MethylStatus"]]
        read_count = len(read_pos_methyl["Read"].unique())

        # Keep only meth positions with counts in at least 15% of all reads in amplicon
        pos_to_keep = meth_pos_counts[meth_pos_counts > read_count*ampl_prev_thr].index
        # posToKeepCount = len(pos_to_keep)
        read_pos_methyl = read_pos_methyl[read_pos_methyl["Pos"].isin(pos_to_keep)]

        # reshape the DF
        read_pos_methyl.set_index(["Read", "Pos"], inplace=True)
        methyl_pattern = read_pos_methyl.unstack()
        methyl_pattern.reset_index(inplace=True)
        methyl_pattern = methyl_pattern.sort_index(axis=1)
        methyl_pattern = methyl_pattern.drop(labels = "Read", axis = 1)
        methyl_pattern.columns = methyl_pattern.columns.droplevel()

        # Fill in NaN with asteriscs (NaN in the case methylation site not on all reads)
        methyl_pattern = methyl_pattern.fillna("*")
        # Collapses identical methylation patterns together and adds  column with the count for each pattern
        collapsed_counted_patterns = methyl_pattern.groupby(
            methyl_pattern.columns.tolist()).size().reset_index().rename(columns={0:'counts'})
        # totalMethPos = methyl_pattern.shape[1]
        # Count the per read methylation states and save in a separate column
        collapsed_counted_patterns["methStatesCount"] = count_states(collapsed_counted_patterns, "+")
        collapsed_counted_patterns["unmethStatesCount"] = count_states(collapsed_counted_patterns, "-")
        collapsed_counted_patterns["notApplCount"] = count_states(collapsed_counted_patterns, "*")

        # Save in table
        p = Path(outpath + "/perSample/")
        p.mkdir( exist_ok=True)
        collapsed_counted_patterns.to_csv("{0}/perSample/{1}_{2}_{3}.{4}.tsv".
                                          format(outpath, sample_id, chrom, snp_coord, allele), sep ="\t", header=True)
        count_methyl_CpGs = pd.pivot_table(collapsed_counted_patterns,
                                    values="counts",
                            index="methStatesCount", aggfunc=np.sum)
        count_methyl_CpGs.rename(columns = {"counts": allele}, inplace=True)

        # Splits the methylation patterns in 3 categories:
        # Mostly methylated (meth >= totalMethPos-methylThr)
        # Mostly unMethylated (meth <= methylThr)
        # Else patrially_meth
        # Count reads in each category
        # Returns a Series
        methylated = sum(collapsed_counted_patterns[
                             collapsed_counted_patterns["methStatesCount"] >= (upper_mCG_thr)]["counts"])
        unmethylated = sum(collapsed_counted_patterns[collapsed_counted_patterns["methStatesCount"] <= methyl_thr]["counts"])
        patrially_meth = read_count - methylated - unmethylated
        methyl_pcnt = methylated / read_count
        unmethyl_pcnt = unmethylated / read_count
        partial_pcnt = patrially_meth / read_count
        count_meth_class = pd.Series(
            [read_count, methylated, methyl_pcnt, unmethylated, unmethyl_pcnt, patrially_meth, partial_pcnt],
            index=["totalReads",
                   "methylated_reads(mCpGs>={})".format(upper_mCG_thr),
                   "methylPcnt",
                   "unmethylated_reads(mCpGs<={})".format(methyl_thr),
                   "unmethylPcnt",
                   "patriallyMeth_reads",
                   "partialPcnt"])
    else:
        count_meth_class = pd.Series(
            [0, 0, 0, 0, 0, 0, 0],
            index=["totalReads",
                   "methylated_reads(mCpGs>={})".format(upper_mCG_thr),
                   "methylPcnt",
                   "unmethylated_reads(mCpGs<={})".format(methyl_thr),
                   "unmethylPcnt",
                   "patriallyMeth_reads",
                   "partialPcnt"])
        # Create an empty DF so the methyl_DFs can be created
        count_methyl_CpGs = pd.DataFrame()
    methyl_DFs = MethylationDFs(count_meth_class, count_methyl_CpGs)

    return methyl_DFs


def count_states(meth_matrix, meth_state):
    """
    Count the occurrence of strings (denoting methylations states) in the table per pattern (row):
    Returns a Series with length equal to matrix rows
    """
    patterns = meth_matrix.drop(labels="counts", axis=1)
    counts = patterns.apply(func = lambda x: sum(x == meth_state), axis = 1)
    return counts
