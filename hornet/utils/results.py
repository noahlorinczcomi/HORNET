# standard imports
import os
import subprocess
import pd as pd
import np as np

# custom imports
from MNet import *  # this is a source file in the HORNET directory


def organizeMetaResults(_dict_):
    dNames = list(_dict_)
    # in a running data frame with _dict_ elements indexed by group ID (ie gene group number)
    df0 = pd.DataFrame.from_dict(_dict_[dNames[0]]["meta"])
    df0["groupID"] = dNames[0]
    for _ in range(1, len(dNames)):
        innerD = _dict_[dNames[_]]["meta"]
        pp = pd.DataFrame.from_dict(innerD)
        pp["groupID"] = dNames[_]
        df0 = pd.concat([df0, pp])
    # df0=df0.reset_index()
    # df0=df0.drop(columns='index', axis=1)
    return df0


def concatUniMRRes(outerDict):
    keys = list(outerDict)
    start = outerDict[keys[0]]
    for II in range(1, len(keys)):
        start = pd.concat((start, outerDict[keys[II]]))
    return start


# organize things to monitor
def organizeThingsMonitored(_dict_):
    df0 = pd.DataFrame.from_dict(_dict_).T
    df0["groupID"] = df0.index
    return df0


# writes out LD matrix to writableDir to be read in by another program (R program I will make) to make plots
def LDLeadGeneSNPs(
    res,
    writableDir=os.getcwd(),
    ldRefDir="/mnt/rstor/SOM_EPBI_XXZ10/njl96/data/1000G/1kg_phase3_EUR_only",
):
    # res: direct output from organizeMetaResults(outerDict)
    nondupres = res.drop_duplicates(subset="LeadSNP", keep="first")
    # note, some genes may share the same lead SNPs ... PLINK will only consider a SNP once.
    outDir1 = (
        writableDir + "/myExtract.txt"
    )  # (for saving SNPs for LD calculations) write out to a directory that must be writable (since the user is there)
    outDir2 = writableDir + "/tempOut"  # (for saving LD matrix)
    res["LeadSNP"].to_csv(
        outDir1, header=False, index=False, sep=" "
    )  # writing out list of SNPs to give PLINK to make LD calculations
    cmd = [
        callPlink(),
        "--r",
        "square",
        "--bfile",
        ldRefDir,
        "--extract",
        outDir1,
        "--out",
        outDir2,
    ]  # the command for PLINK
    _ = subprocess.call(
        cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
    )  # executing the command in the terminal
    ldmat = np.loadtxt(
        writableDir + "/tempOut.ld", dtype="f"
    )  # reading LD matrix back in; matrix format but with duplicates occupying only one row/column
    ### maybe for the plotting this does not matter actually ...
    # ldmat is for unique SNPs only. write out correlation matrix for nondupres, which is all SNPs but none repeated >once, and the group start indices
    # if writing out for another program to read in
    os.chdir(writableDir)
    np.savetxt("ldmat.txt", ldmat)
    res.to_csv("res.txt", index=False)
    nondupres.to_csv("nondupres.txt", index=False)
    # note what this function does in the end - it doesnt return anything; it writes out a file


def prepRes(res):
    res = res[
        [
            "Gene",
            "geneBP",
            "GeneSelected",
            "MRBEEPosSelEst",
            "MRBEEPosSelSE",
            "LocusR2",
            "geneZ",
            "geneSNP",
            "snpBP",
            "conditionalF",
            "groupID",
            "CHR",
            "infl",
            "h2CisEstNaive",
            "nClumps",
            "SingleSNPEst",
            "SingleSNP_P",
            "IVW_MVMR_Est",
            "IVW_MVMR_SE",
            "IVW_UNIMR_Est",
            "IVW_UNIMR_SE",
            "UniMRBEEEsts",
            "UniMRBEESEs",
            "MRBEEEsts",
            "MRBEESEs",
            "SMRm",
            "SMRpleioJointPs",
            "SMRHPmeansP",
            "OVBiasTestP",
        ]
    ]
    res = res.rename(
        columns={
            "geneZ": "LeadEQTLZstatistic",
            "geneSNP": "leadEQTL",
            "snpBP": "leadEQTLBP",
            "groupID": "CHRspecificGroupID",
            "CHR": "Chromosome",
            "infl": "nullInflation",
            "IVW_UNIMR_Est": "IVW_UVMR_Est",
            "IVW_UNIMR_SE": "IVW_UVMR_SE",
            "UniMRBEEEsts": "MRBEE_UVMR_Est",
            "UniMRBEESEs": "MRBEE_UVMR_SE",
            "MRBEEEsts": "MRBEE_MVMR_Est",
            "MRBEESEs": "MRBEE_MVMR_SE",
            "MRBEEPosSelEst": "MRBEEPostSelec_MVMR_Est",
            "MRBEEPosSelSE": "MRBEEPostSelec_MVMR_SE",
            "SMRm": "SMRNumIVs",
            "SMRHPmeansP": "SMRUnbalancedHP_P",
            "OVBiasTestP": "SMROVBias_P",
        }
    )
    return res
