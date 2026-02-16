# standard imports
from scipy.stats import norm
import platform
import subprocess
import pandas
import numpy as np

# custom imports
from MNet import *  # this is a source file in the HORNET directory


def callPlink():
    sys = platform.system()
    if sys == "Linux":
        call_ = "./plinkdir/linux/plink"
    elif sys == "Darwin":
        call_ = "./plinkdir/mac/plink"
    elif sys == "Windows":
        call_ = "plinkdir/windows/plink.exe"
    return call_


def findOutcomeSignals(
    dataPheno,
    ldRefDir,
    writableDir,
    outcomeClumpingKBWindow,
    outcomeClumpingPthreshold,
    outcomeClumpingR2,
):
    # outcomeClumpingKBWindow, outcomeClumpingPthreshold, outcomeClumpingR2:
    #  - parameters to give PLINK to find clumps in the outcome GWAS (whole genome, not CHR-specific)
    dataPheno_ = dataPheno.copy()[
        abs(dataPheno["phenoZ"]) > norm.ppf(1 - outcomeClumpingPthreshold)
    ]  # subset to speed things up
    dataPheno_["P"] = np.array(
        [
            norm.cdf(-abs(dataPheno_["phenoZ"].values[x]))
            for x in range(0, dataPheno_.shape[0])
        ]
    )
    dataPheno_[["phenoSNP", "P"]].rename(columns={"phenoSNP": "SNP"}).to_csv(
        writableDir + "/outcomePsout.txt", index=False, sep="\t"
    )
    cmd = [
        callPlink(),
        "--bfile",
        ldRefDir,
        "--clump",
        writableDir + "/outcomePsout.txt",
        "--clump-kb",
        str(outcomeClumpingKBWindow),
        "--clump-p1",
        str(outcomeClumpingPthreshold),
        "--clump-p2",
        str(outcomeClumpingPthreshold),
        "--clump-r2",
        str(outcomeClumpingR2),
        "--out",
        writableDir + "/outcomeplinkout",
    ]
    out = subprocess.call(
        cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
    )  # executing the command in the terminal
    outcomeClumps = pandas.read_fwf(
        writableDir + "/outcomeplinkout.clumped", sep="\t"
    )  # reading LD matrix back in
    outcomeClumps = (
        outcomeClumps["SNP"].dropna().values
    )  # for some reason some NAs get appended to the end
    dataPheno["isOutcomeClump"] = dataPheno["phenoSNP"].isin(outcomeClumps)
    return dataPheno
