# standard imports
from scipy.stats import norm
import subprocess
import pandas
import numpy as np

# custom imports
from MNet import *  # this is a source file in the HORNET directory


def UVMRworkhorse(
    merged,
    writableDir,
    ldRefDir,
    LDWindowMb=1,
    LDMaxInWindow=0.5,
    eQTLInOtherWindowP=5e-5,
    IVP=0.025,
    minNumIVs=10,
    LDWithinIVSetMax=0.9,
    restrictToOutcomeSignals=True,
    outcomeSignalsWithinXMb=1,
):
    # LDWindowMb: remove IVs that have LD>LDMaxInWindow with other eQTLs within LDWindowMb
    # LDMaxInWindow: IVs with LD>LDMaxInWindow with other eQTLs within LDWindowMb are removed
    # eQTLInOtherWindowP: from two above, only consider SNPs in nearby windows to be eQTLs for other genes if P<eQTLInOtherWindowP
    # IVP: P-value threshold for defining IVs in UniMR
    # minNumIVs: self-explanatory
    genes = merged["Gene"].unique()
    outcomeClumpBPs = merged[merged["isOutcomeClump"] == True]["snpBP"].values
    mergedsig = merged[abs(merged["geneZ"]) > norm.ppf(1 - eQTLInOtherWindowP)]
    thingsToMonitor = {}
    outerDict = {}
    for ogene in range(0, len(genes)):
        mcut = merged.copy()[merged["Gene"] == genes[ogene]]
        mcut = mcut[abs(mcut["geneZ"]) > norm.ppf(1 - IVP)]
        if mcut.shape[0] < minNumIVs:  # if not enough eQTLs, next gene
            continue
        # if user only wants to do the analysis near outcome signals, make sure there are some here
        if restrictToOutcomeSignals == True:
            diffs = [
                min(abs((mcut["snpBP"].values)[x] - outcomeClumpBPs))
                for x in range(0, mcut.shape[0])
            ]
            if any(np.array(diffs) < (1e6 * outcomeSignalsWithinXMb)) == False:
                continue
        # find LD between these SNPs
        mcut = mcut.drop_duplicates(
            subset=["geneSNP"], keep=False
        )  # drop duplicate SNPs (cannot appear in file given to PLINK)
        mcut = mcut.sort_values(by="snpBP", ascending=True)
        if mcut.shape[0] < minNumIVs:  # if no non-duplicated eQTLs, next gene
            continue
        # make sure not in LD with eQTLs with any genes within LDWindowMb
        minBP, maxBP = min(mcut["snpBP"]), max(mcut["snpBP"])
        mask1 = mergedsig[
            (mergedsig["snpBP"] < minBP)
            & (mergedsig["snpBP"] > (minBP - LDWindowMb * 1e6))
            & (mergedsig["Gene"] != genes[ogene])
        ]
        mask2 = mergedsig[
            (mergedsig["snpBP"] > maxBP)
            & (mergedsig["snpBP"] < (maxBP + LDWindowMb * 1e6))
            & (mergedsig["Gene"] != genes[ogene])
        ]
        mask1 = (
            mask1[["geneSNP", "snpBP"]]
            .drop_duplicates(keep="first")
            .sort_values(by="snpBP", ascending=True)
        )
        mask2 = (
            mask2[["geneSNP", "snpBP"]]
            .drop_duplicates(keep="first")
            .sort_values(by="snpBP", ascending=True)
        )
        # write out big list in this order: mask1, reshaped, mask2
        big_list = pandas.concat([mask1["geneSNP"], mcut["geneSNP"], mask2["geneSNP"]])
        outDir1 = writableDir + "/myExtract.txt"
        outDir2 = writableDir + "/ld"
        big_list.to_csv(
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
        out = subprocess.call(
            cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT
        )  # executing the command in the terminal
        ldmat = np.loadtxt(outDir2 + ".ld", dtype="f")  # reading LD matrix back in
        # now need to know if any SNPs in reshaped are in LD>ldOtherLociR2(=0.1 be default) with SNPs outside of reshaped
        maskreshaped = np.ones((ldmat.shape[0],), dtype="bool")
        maskreshaped[: mask1.shape[0]] = False
        maskreshaped[(mask1.shape[0] + mcut.shape[0]) :] = False
        inds1 = ldmat[mask1.shape[0] : (mask1.shape[0] + mcut.shape[0]), :][
            :, : mask1.shape[0]
        ]
        inds1 = np.nan_to_num(inds1, nan=0)
        inds2 = ldmat[mask1.shape[0] : (mask1.shape[0] + mcut.shape[0]), :][
            :, (mask1.shape[0] + mcut.shape[0]) :
        ]
        inds2 = np.nan_to_num(inds2, nan=0)
        ldmask1 = (((inds1**2) > (LDMaxInWindow**2)).sum(axis=1)) > 0
        ldmask2 = (((inds2**2) > (LDMaxInWindow**2)).sum(axis=1)) > 0
        ldmask = (ldmask1 + ldmask2) == 0
        if sum(ldmask) < minNumIVs:  # if no SNPs that are valid IVs for univariable MR
            continue
        mcut = mcut[ldmask]
        ldmat = ldmat[maskreshaped, :][:, maskreshaped]  # subset to only SNPs in mcut
        ldmat = ldmat[ldmask, :][
            :, ldmask
        ]  # subset again to only SNPs in reshaped uncorrelated with SNPs in surrounding window
        # some LD correlations are totally missing for a SNP (rare, but possible). remove these. (possible why? observed AF=0 or 1 so no variation?)
        mask = (
            (np.isnan(ldmat).sum(axis=0)) == ldmat.shape[0]
        ) == False  # create mask: True if non-missing LD; False if all missing LD (I've never seen partial missing)
        mcut = mcut[mask]
        del mask1, mask2, ldmask, inds1, inds2, ldmask1, ldmask2
        someDropped, snpsKept = LDShrinkUpper(
            ldmat, LDWithinIVSetMax, mcut["geneSNP"].values
        )
        ld0 = someDropped.copy()  # find those not in LD too highly
        mcut = mcut[mcut["geneSNP"].isin(snpsKept)]
        # make sure LD between remaining SNPs isn't too high to cause inflation
        if ld0.shape[0] < minNumIVs:
            continue
        # perform MR (bias-correction for SigmaUV assumed to be 0 since very unlikely gene expression samples overlap with any phenotype samples)
        # first step is to transform by sqrt of inverse of LD matrix
        bxog, byog = mcut["geneZ"].values, mcut["phenoZ"].values
        m = bxog.shape[0]
        bxog, byog = bxog.reshape((m, 1)), byog.reshape((m, 1))
        # but first make LD matrix positive definite
        ld0, alpha, deco = regularizeLD(
            ld0, int(ld0.shape[0] / 4), ld0.shape[0] + 1, 0.1, epsilon=1e-4, toCorr=True
        )  # int() rounds down
        ld0inv = deco[1] @ np.diag(1 / deco[0]) @ deco[1].T
        ld0invsqrt = deco[1] @ np.diag(1 / deco[0] ** 0.5) @ deco[1].T
        bx = ld0invsqrt @ bxog
        by = ld0invsqrt @ byog
        ld0 = np.eye(m)
        ld0inv = ld0.copy()
        ld0invsqrt = ld0.copy()
        bx_ = np.column_stack(([1] * m, bx))
        UU = np.ones((1, 1))  # under transformation, variance of measurement error is 1
        UV = np.zeros((1, 1))
        VV = np.ones((1, 1))
        UU_ = np.zeros((2, 2))
        UU_[1, 1] = 1
        UV_ = np.zeros((2, 1))
        # MRBEE
        est0, V0, Outliers, kR, k = imrbee(
            bx,
            by,
            UU,
            UV,
            VV,
            ld0inv,
            ld0,
            PleioPThreshold=0.05 / bx.shape[0] ** 0.5,
            boot=False,
            niter=1000,
            initial="bee",
            max_iter=15,
        )
        pleiop = SpleioP(bx_, by, UU_, UV_, VV, est0, V0)
        # IVW (MR-Egger actually)
        bxx = np.linalg.inv(bx_.T @ bx_)
        ivwest = bxx @ bx_.T @ by
        s2 = sum((by.squeeze() - bx_ @ ivwest.squeeze()) ** 2) / (m - 2)
        ivwvar = s2 * bxx
        # huber loss
        huberE, huberV, huberW = MRBEEHuber(
            bx,
            by,
            ld0,
            UU,
            UV,
            gamma=0.5,
            initial="bee",
            eps=0.001,
            max_iter=15,
            boot=False,
        )  # keep gamma small (as gamma-> huber->quantile)
        # end
        paraDict = {
            "Gene": genes[ogene],
            "GeneBP": int(mcut["geneBP"].unique()),
            "MRBEEEst": float(est0[1]),
            "MRBEESE": float(V0[1, 1] ** 0.5),
            "HuberEst": float(huberE[1]),
            "HuberSE": float(huberV[1, 1] ** 0.5),
            "IVWEst": float(ivwest[1]),
            "IVWSE": float(ivwvar[1, 1] ** 0.5),
            "nIVs": m,
        }
        outerDict[genes[ogene]] = pandas.DataFrame(paraDict, index=[0])
    return outerDict
