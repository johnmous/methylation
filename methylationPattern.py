## Author: I. Moustakas
## Title: Get the methylation patterns and count them
## Usage: methylationPattern.py bismarkCpG outputDir sampleName

import click
import os
import numpy as np
import pandas as pd
import re


##
## @click.command()
## @click.option('--cpgfile', required = True, help = "CpG file from bismark (CpG_OB_*)")
## @click.option('--ampltable', required = True, help = "Tab separated file with amplicon locations") ## TODO: specify table format
## @click.option('--outpath', required = True, help='Path to place the output')


def methylPatterns(methylation, outpath):
    """
    Separate reads in three categories: Methylated , unmethylated, partially methylated.
    Count each category.
    :param methylation: pandas df with columns "Read", "MethylStatus", "Chr", "Pos", "Zz"
    :param ampltable:
    :param outpath:
    :return: A Series with counts of the three methylation categories per amplicon
    """

    ## if records for amplicon, proceeed
    if len(methylation.index)>0:
        ## Get all methylation positions in the amplicon
        methPosCounts = methylation["Pos"].value_counts()
        readPosMethyl = methylation[["Read", "Pos", "MethylStatus"]]
        readCount = len(readPosMethyl["Read"].unique())

        ## Keep only meth posistions with counts in at least 1% of all reads in amplicon
        posToKeep = methPosCounts[methPosCounts > readCount*0.01].index
        posToKeepCount = len(posToKeep)
        readPosMethyl = readPosMethyl[readPosMethyl["Pos"].isin(posToKeep)]

        ## reshape the DF
        readPosMethyl.set_index(["Read", "Pos"], inplace=True)
        methylPattern = readPosMethyl.unstack()
        methylPattern.reset_index(inplace=True)
        methylPattern = methylPattern.drop(labels = "Read", axis = 1)
        methylPattern.columns = methylPattern.columns.droplevel()

        ## Fill in NaN with asteriscs (NaN in the case methylation site not on all reads)
        ## Count the methylation patterns
        methylPattern = methylPattern.fillna("*")
        collCountPattern = methylPattern.groupby(methylPattern.columns.tolist()).size().reset_index().rename(columns={0:'counts'})
        totalMethPos = methylPattern.shape[1]

        ## Count the per pattern methylation states
        methStates = countStates(collCountPattern, "+")
        unmethStates = countStates(collCountPattern, "-")
        notAppl = countStates(collCountPattern, "*")
        collCountPattern["methStatesCount"] = methStates
        collCountPattern["unmethStatesCount"] = unmethStates
        collCountPattern["notApplCount"] = notAppl
        collCountPattern.to_csv(outpath + "/test.tsv", sep ="\t", header=True)
        countMethClass = countPatterns(collCountPattern, totalMethPos, readCount, 3, posToKeepCount)
        return countMethClass

def countStates(methMatrix, methState):
    """
    Count the number of methylation states in the table per pattern (row):
    methylated (+),
    unmethylated (-),
    not applicable (*)
    Returns a Series with length equal to matrix rows
    """
    patterns = methMatrix.drop(labels="counts", axis=1)
    counts = patterns.apply(func = lambda x: sum(x == methState), axis = 1)
    return counts

def countPatterns(methMatrix, totalMethPos, readCount, stateThr, posToKeepCount):
    """
    Split the methylation patterns in 3 categories:
    Mostly methylated (meth > stateThr)
    Mostly unMethylated (meth < stateThr)
    Else patriallyMeth
    Count reads in each category
    Returns a Series
    """
    methylated = sum(methMatrix[methMatrix["methStatesCount"]>totalMethPos-stateThr]["counts"])
    unmethylated = sum(methMatrix[methMatrix["unmethStatesCount"]>stateThr]["counts"])
    patriallyMeth = readCount - methylated - unmethylated
    methylPcnt = methylated/readCount
    unmethylPcnt = unmethylated/readCount
    partialPcnt = patriallyMeth/readCount
    countMethClass = pd.Series([readCount, methylated, methylPcnt, unmethylated, unmethylPcnt, patriallyMeth, partialPcnt],
                             index = ["totalReads",
                                      "methylated_reads_{}-{}mGC".format(posToKeepCount, posToKeepCount-stateThr),
                                      "methylPcnt",
                                      "unmethylated_reads_{}-{}".format(posToKeepCount, posToKeepCount-stateThr),
                                      "unmethylPcnt",
                                      "patriallyMeth_reads",
                                       "partialPcnt"])
    return countMethClass

def sampleName(file):
    """
    Get sample name out of the bismakr file name.
    Expects full path of a file, as named by bismark
    """
    str_search = re.search('.+/CpG_OB_(.+)_bismark.+', file)
    sampleName = str_search.group(1)
    return sampleName



if __name__ == "__main__":
    methylPatterns()
