import pandas as pd
from pathlib import Path


def methyl_patterns(methyl_extr, outpath, methyl_thr, number_CGs, sample_id, allele, chrom, snp_coord):
    """
    Separate reads in three categories: Methylated , unmethylated, partially methylated according to the number
    of CpGs that are methyalted
    :param methyl_extr: pandas df with columns "Read", "MethylStatus", "Chr", "Pos", "Zz". This is the output of bismark_methylation_extractor
    :param outpath: path to save the table with the methylation patterns
    :param: number_CGs: Number of CGs in the amplicon
    :param methyl_thr: Integer to subtract from number_CGs to set threshold from mostly methylated
    :return: A pandas.Series with read counts and percentages for the three methylation categories
    """

    # if records for amplicon
    if len(methyl_extr.index)>0:

        # Get all methylation positions in the amplicon
        meth_pos_counts = methyl_extr["Pos"].value_counts()
        read_pos_methyl = methyl_extr[["Read", "Pos", "MethylStatus"]]
        read_count = len(read_pos_methyl["Read"].unique())

        # Keep only meth posistions with counts in at least 1% of all reads in amplicon
        pos_to_keep = meth_pos_counts[meth_pos_counts > read_count*0.01].index
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

        # Splits the methylation patterns in 3 categories:
        # Mostly methylated (meth >= totalMethPos-methylThr)
        # Mostly unMethylated (meth <= methylThr)
        # Else patrially_meth
        # Count reads in each category
        # Returns a Series
        methylated = sum(collapsed_counted_patterns[collapsed_counted_patterns["methStatesCount"] >= (number_CGs - methyl_thr)]["counts"])
        unmethylated = sum(collapsed_counted_patterns[collapsed_counted_patterns["methStatesCount"] <= methyl_thr]["counts"])
        patrially_meth = read_count - methylated - unmethylated
        methyl_pcnt = methylated / read_count
        unmethyl_pcnt = unmethylated / read_count
        partial_pcnt = patrially_meth / read_count
        count_meth_class = pd.Series(
            [read_count, methylated, methyl_pcnt, unmethylated, unmethyl_pcnt, patrially_meth, partial_pcnt],
            index=["totalReads",
                   "methylated_reads(mGCs>={})".format(number_CGs - methyl_thr),
                   "methylPcnt",
                   "unmethylated_reads(mGCs<={})".format(methyl_thr),
                   "unmethylPcnt",
                   "patriallyMeth_reads",
                   "partialPcnt"])
    else:
        count_meth_class = pd.Series(
            [0, 0, 0, 0, 0, 0, 0],
            index=["totalReads",
                   "methylated_reads(mGCs>={})".format(number_CGs - methyl_thr),
                   "methylPcnt",
                   "unmethylated_reads(mGCs<={})".format(methyl_thr),
                   "unmethylPcnt",
                   "patriallyMeth_reads",
                   "partialPcnt"])

    return count_meth_class


def count_states(meth_matrix, meth_state):
    """
    Count the occurence of strings (denoting methylations states) in the table per pattern (row):
    Returns a Series with length equal to matrix rows
    """
    patterns = meth_matrix.drop(labels="counts", axis=1)
    counts = patterns.apply(func = lambda x: sum(x == meth_state), axis = 1)
    return counts
