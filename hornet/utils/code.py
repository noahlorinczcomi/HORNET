# standard imports
import os
import platform
import pandas
import numpy as np

# custom imports
from MNet import *  # this is a source file in the HORNET directory


def makeI(k):  # function to define a kxk identity matrix
    I = np.zeros((k, k))
    np.fill_diagonal(I, 1)
    return I


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def unpackDictOfLists(ddict):
    els = []
    for i in range(0, len(ddict)):
        els.append(list(ddict[list(ddict)[i]]))
    return flatten_list(els)


def quickAppend(nestedlist):
    out = [nestedlist[x].append(float("nan")) for x in range(0, len(nestedlist))]
    return out


def c(A):  # function to vectorize (going down columns) a matrix A
    return t(np.asarray(A)).ravel()


def whichMax(ll):
    return list(np.where(np.array(ll) == max(ll)))[0][0]


def reshapeCorrelDf(merged, theseGenes, geneLabel, joinType):
    dd = {}
    newColNames = []
    for jj in range(0, len(theseGenes)):
        dfGene = merged[merged[geneLabel] == theseGenes[jj]]
        dfGene = dfGene[["geneSNP", "geneZ"]]
        # rename geneZ to be gene-specific
        cn = theseGenes[jj] + "_Z"
        newColNames.append(cn)
        dfGene = dfGene.rename(columns={"geneZ": cn})
        dd[theseGenes[jj]] = dfGene
        if jj > 0:
            prevkey = list(dd)[jj - 1]
            prevdat = dd[prevkey]
            if joinType == "inner":
                binding = pandas.merge(
                    prevdat, dfGene, how="inner", left_on="geneSNP", right_on="geneSNP"
                )
            else:
                binding = pandas.merge(
                    prevdat, dfGene, how="outer", left_on="geneSNP", right_on="geneSNP"
                )
            dd[theseGenes[jj]] = binding
    # now add phenotype to this, first column
    toBindAgain = merged[["phenoSNP", "phenoZ"]].drop_duplicates(["phenoSNP", "phenoZ"])
    binding0 = pandas.merge(
        toBindAgain, binding, left_on="phenoSNP", right_on="geneSNP"
    )
    mergedcut = merged.drop_duplicates(subset=["geneSNP"], keep="first")
    binding0 = pandas.merge(
        mergedcut.drop(columns=["geneZ", "phenoZ"]),
        binding0,
        left_on="geneSNP",
        right_on="geneSNP",
    )
    return binding0


def flatten_list(_2d_list):  # flatten list of lists into a single list
    flat_list = []
    # Iterate through the outer list
    for element in _2d_list:
        if type(element) is list:
            # If the element is of type list, iterate through the sublist
            for item in element:
                flat_list.append(item)
        else:
            flat_list.append(element)
    return flat_list


def subR(R, removeInds):
    a1 = np.delete(R, removeInds, 1)
    a2 = np.delete(a1, removeInds, 0)
    return a2


def makeMask(n, inds, starting=True):
    mask = (
        np.ones((n,), dtype="bool")
        if starting == True
        else np.zeros((n,), dtype="bool")
    )
    mask[inds] = starting == False
    return mask


def indexFinder(toFind, allValues):
    gg = [
        allValues.index(toFind[x])
        for x in range(0, len(toFind))
        if toFind[x] in allValues
    ]
    return gg


def columnFinder(columnNames):
    columnNames = [x.lower() for x in columnNames]
    # columnNames must be a list of column names in either the gene expression or outcome data set
    # goal is to auto-detect what columns in the data mean
    # columns to detect: SNP, CHR, BP, A1, A2, Z-statistics
    # potential names
    snpNames = [
        "rsid",
        "snp",
        "markername",
    ]  # need to make sure not in CHR:BP:A1:A2 format
    bpNames = [
        "pos",
        "bp",
        "bp_hg19",
        "bphg19",
        "pos_hg19",
        "poshg19",
        "snppos",
        "snpbp",
    ]
    chrNames = ["chr", "chrom", "chr_h19", "chrhg19", "snpchr", "snpchrom"]
    a1Names = [
        "a1",
        "allele1",
        "effectallele",
        "effect_allele",
        "assessed_allele",
        "assessedallele",
        "ref",
    ]
    a2Names = [
        "a2",
        "allele2",
        "noneffectallele",
        "noneffect_allele",
        "nonassessed_allele",
        "assessedallele",
        "other_allele",
        "otherallele",
        "ref",
    ]
    zNames = [
        "z",
        "zstat",
        "zstatistic",
        "z_stat",
        "z_statistic",
        "t",
        "tstat",
        "t_stat",
        "t_statistic",
        "teststatistic",
        "test_statistic",
        "zscore",
        "z_score",
        "tscore",
        "t_score",
    ]
    # find their indices
    snpInd = indexFinder(snpNames, columnNames)
    bpInd = indexFinder(bpNames, columnNames)
    chrInd = indexFinder(chrNames, columnNames)
    a1Ind = indexFinder(a1Names, columnNames)
    a2Ind = indexFinder(a2Names, columnNames)
    zInd = indexFinder(zNames, columnNames)
    # double-check all were found
    keys = ["SNP", "BP", "CHR", "A1", "A2", "Z"]
    vals = [snpInd, bpInd, chrInd, a1Ind, a2Ind, zInd]
    if len(flatten_list(vals)) < len(keys):
        raise ValueError(
            "Could not find column name for " + keys[vals.index([])] + " column"
        )
    # make dictionary
    dict = {
        "snpInd": snpInd[0],
        "bpInd": bpInd[0],
        "chrInd": chrInd[0],
        "a1Ind": a1Ind[0],
        "a2Ind": a2Ind[0],
        "zInd": zInd[0],
    }
    return dict


def bimin(mat):
    return np.where(mat == mat.min())


def callDelete():
    sys = platform.system()
    if (sys == "Linux") | (sys == "Darwin"):
        call_ = "rm "
    elif sys == "Windows":
        call_ = "del "
    return call_


def fileChecker(fp, ftype):
    if os.path.exists(fp) == False:
        raise ValueError(
            "The {} does not exist. Are you sure you entered it correctly?".format(
                ftype
            )
        )
