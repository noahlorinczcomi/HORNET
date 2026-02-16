# standard imports
from scipy.stats import norm
import scipy
import subprocess
from scipy import stats
import numpy as np
import math

# custom imports
from MNet import *  # this is a source file in the HORNET directory


def fPos(x):
    cc = x < 0
    x[cc] = 0
    return x


def mad(a):
    b = 1.483 * np.median(abs(a - np.median(a)))
    return b


def chiTestInternal(v, R):
    v = v.reshape((v.shape[0], 1))
    cc_ = v.T @ np.linalg.inv(R) @ v
    cc_ = cc_[0]
    cp_ = 1 - stats.chi2.cdf(cc_, v.shape[0])
    return cp_


def IVChiTest(Zmat, R):
    Z0 = np.nan_to_num(Zmat, nan=0)
    chiPs = []
    for _ in range(0, Z0.shape[0]):
        ta = chiTestInternal(Z0[_, :], R)
        chiPs.append(ta)
    chiPs = np.array(chiPs)
    return chiPs


def CD_BIC(
    x,
    y,
    UU,
    UV,
    VV,
    Rinv,
    lambda_min,
    lambda_max,
    alpha_,
    nLambda=100,
    eps=0.00001,
    max_iter=100,
    norm_stoppage=100,
    pen="scad",
):
    # x: bx0, y: by, UU: SigmaUU, UV: SigmaUV, VV:SigmaVV, lambda_min/_max: minimum/maximum lambda, alpha_: alpha
    lambdas = np.linspace(lambda_min, lambda_max, nLambda)
    bic = []
    for ll in range(0, lambdas.shape[0]):
        est_ll = CD(
            x,
            y,
            UU,
            UV,
            VV,
            Rinv,
            lambdas[ll],
            alpha_,
            eps,
            max_iter,
            norm_stoppage,
            pen,
        )["est"]
        mse = sum((y - x @ est_ll) ** 2)[0]
        # df is trace of hat matrix, but that can be computationally costly since it is one large matrix inversion per lambda ... stick with #nonzero
        # Hat=x@np.linalg.inv(x.T@Rinv@x-m*UU)@x.T
        # HatRidge=x@np.linalg.inv(x.T@Rinv@x-m*UU+lambdas[ll]*makeI(p))@x.T
        # dfRidge=np.trace(HatRidge)
        # df=np.trace(Hat) # not even clear what df are for MRBEE
        effDF = len([x for x in range(0, est_ll.shape[0]) if est_ll[x][0] != 0])
        bic.append(np.log(mse) / m + effDF * np.log(m) / m)
    dd = {"BICs": bic, "lambda_min": lambdas[bic.index(min(bic))]}
    return dd


def nuclearNorm(A):
    return sum(np.linalg.svd(A)[1])


def poet(A, k):
    u, s, vh = np.linalg.svd(A)
    u = u[:, :k]
    return u


def ldscore(ld0):
    ell = []
    mask = np.zeros((ld0.shape[0],), dtype="bool")
    for _o in range(0, ld0.shape[0]):
        mask_o = mask.copy()
        mask_o[_o] = True
        _v = ld0[mask_o, :][:, (mask_o == False)].reshape((ld0.shape[0] - 1,))
        ell.append(sum(_v**2))
    ell = np.array(ell)
    return ell


def wmean(x, w):
    return sum(x * w) / sum(w)


def TwoMixtureEM(data, mu01, mu02, sig01, sig02, eps=0.001, k=100):
    p01 = norm.pdf(data, mu01, sig01)
    p02 = norm.pdf(data, mu02, sig02)
    error = 1
    k_ = 0
    es = []
    gg = np.ones((data.shape[0]), dtype="int")
    while (error > eps) & (k_ < k):
        k_ = k_ + 1
        mask = p01 < p02  # True if predicted group 2, False if predicted group 1
        mu1est = wmean(
            data[mask == False], p01[mask == False]
        )  # can change weights if you want
        mu2est = wmean(data[mask], p02[mask])
        sig1est = np.std(x[mask == False])
        sig2est = np.std(x[mask])
        p01 = norm.pdf(data, mu1est, sig1est)
        p02 = norm.pdf(data, mu2est, sig2est)
        error = (mu1est - mu01 + mu2est - mu02 + sig1est - sig01 + sig2est - sig02) ** 2
        mu01 = mu1est
        mu02 = mu2est
        sig10 = sig1est
        sig20 = sig2est  # should all be floats so no need to use copies
        es.append(error)
        if k_ > 1:
            if (
                abs(es[-2] - es[-1]) < 0.0001
            ):  # another stopping rule - if error did not change much at this iteration
                error = 0
    gg[mask] = 2
    gg = gg - 1
    return gg


def trace(A):
    out = np.sum(np.diag(A))
    return out


def adjustInflation(res, infl):
    # want to adjust (Est/SE), (IVW_MVMR_Est/IVW_MVMR_SE), (IVW_UNIMR_Est/IVW_UNIMR_SE), (UniMRBEEEsts/UniMRBEESEs),
    # (MRBEEEsts/MRBEESEs), (MRBEERidgeEst/MRBEERidgeSE), (MRBEEPosSelEst/MRBEEPosSelSE)
    res2 = res.copy()
    # res2['MRJonesSE']=res2['MRJonesSE']*(infl**0.5) # I removed MR-Jones's SE so don't need this anymore
    res2["IVW_MVMR_SE"] = res2["IVW_MVMR_SE"] * (infl**0.5)
    res2["IVW_UVMR_SE"] = res2["IVW_UVMR_SE"] * (infl**0.5)
    res2["MRBEE_UVMR_SE"] = res2["MRBEE_UVMR_SE"] * (infl**0.5)
    res2["MRBEE_MVMR_SE"] = res2["MRBEE_MVMR_SE"] * (infl**0.5)
    # res2['MRBEERidge_MVMR_SE']=res2['MRBEERidge_MVMR_SE']*(infl**0.5)
    res2["MRBEEPostSelec_MVMR_SE"] = res2["MRBEEPostSelec_MVMR_SE"] * (infl**0.5)
    # res2['MRBEEHuber_MVMR_SE']=res2['MRBEEHuber_MVMR_SE']*(infl**0.5)
    return res2


def bootVarWLS(x, y, w, UU, k=1000):
    m = len(y)
    p = x.shape[1]
    ests = np.zeros((k, p))
    for ii in range(0, k):
        ind = np.random.randint(0, m, m)
        xi = x[ind]
        yi = y[ind]
        wi = w[ind]
        D = np.diag(wi.squeeze())
        ei = np.linalg.inv(xi.T @ D @ xi - sum(wi) * UU) @ xi.T @ D @ yi
        ests[ii, :] = ei.squeeze()
    return np.cov(ests, rowvar=False)


def ispd(x):
    return np.all(np.linalg.eigvals(x) > 0)


def calcInflationFactor(
    merged,
    ldRefDir,
    writableDir,
    mbWindow=5,
    ldUpperLimit=0.5,
    outcomeP=0.01,
    nGenesToTry=50,
    geneDistMb=1,
    minM=25,
):
    # 1) find loci in which there is no outcome signal
    # 2) select a gene with no eQTLs
    # 3) calculate LD for SNPs available for gene
    # 4) prune LD to only |r|<ldUpperLimit (default=0.5)
    # 5) calculate test statistic and store
    # 6) repeat for many genes that are at least geneDistMb Mb away from each other (default is 5)
    # 7) store inflation
    outDir1 = (
        writableDir + "/myExtract.txt"
    )  # (for saving SNPs for LD calculations) write out to a directory that must be writable (since the user is there)
    outDir2 = writableDir + "/tempOut"  # (for saving LD matrix)
    # 1)
    merged_ = merged.copy().reset_index()
    clumpBPs = merged_[merged_["isOutcomeClump"] == True]["snpBP"].unique()
    mask = np.ones((merged_.shape[0]), dtype="bool")
    for _ in range(0, len(clumpBPs)):
        m_ = merged[
            abs(merged["snpBP"] - clumpBPs[_]) < (mbWindow * 1000000)
        ]  # the 1Mb region around clumpBPs
        mask[m_.index] = False
    merged_ = merged_.sort_values(by="snpBP", ascending=True)
    merged_ = merged_[mask]
    # 2)
    # first finding genes that are at least geneDistMb Mb away from each other
    genedf = (
        merged_[["Gene", "geneBP"]]
        .drop_duplicates()
        .sort_values(by="geneBP", ascending=True)
    )
    if genedf.shape[0] == 0:
        return 1, 0, 1, [], 0, float("nan")
    mask = np.zeros((genedf.shape[0],), dtype="bool")
    mask[0] = True
    for _ in range(1, len(mask)):  # starting with second gene
        tind = np.where(mask == True)[0][-1]
        toadd = (
            True
            if abs(genedf.iloc[_, 1] - genedf.iloc[tind, 1]) > (geneDistMb * 1e6)
            else False
        )
        mask[_] = toadd
        if sum(mask) == nGenesToTry:  # if I already have enough, don't need to continue
            break

    ggs = genedf[mask]["Gene"].unique()
    cstats = list()
    ccs = list()
    ccps = list()
    counter = 0
    usedGenes = list()
    ms = list()
    for _ in range(0, len(ggs)):
        dat = merged_[merged_["Gene"] == ggs[_]]
        dat = dat[abs(dat["phenoZ"]) < norm.ppf(1 - outcomeP, 0, 1)]
        # dat=dat[abs(dat['geneZ'])<norm.ppf(1-outcomeP,0,1)] # also make sure no eQTLs to rule out horizontal pleiotropy causing inflation
        if dat.shape[0] < (minM):
            continue
        # 3) only use every kth SNP to speed things up and ensure I don't need to do as much removing bc of high LD later.
        # I want approximately 100 SNPs before LD pruning so the k in kth is dynamically chosen based on how much data is originally available
        roundFactor = math.floor(dat.shape[0] / 100)
        if roundFactor == 0:
            roundFactor = 1
        dat = dat.iloc[::roundFactor, :].drop_duplicates("geneSNP")
        dat["geneSNP"].to_csv(
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
        ldmat = np.genfromtxt(outDir2 + ".ld", dtype="f")  # reading LD matrix back in
        mask = (
            np.isnan(ldmat).sum(axis=0) == ldmat.shape[0]
        ) == False  # if any missing (bc MAF==0), will have all nan's for that SNP
        ldmat = ldmat[mask, :][:, mask]
        dat = dat[np.array(mask).astype(bool)]
        # 4) prune
        prunedLD, snpsKept = LDShrinkUpper(ldmat, ldUpperLimit, dat["geneSNP"].values)
        dat = dat[dat["geneSNP"].isin(snpsKept)]
        # 5) MR + test stats
        m = dat.shape[0]
        bx = np.column_stack(([1] * m, dat["geneZ"].values))
        by = dat["phenoZ"].values
        res = scipy.stats.spearmanr(dat[["geneZ", "phenoZ"]])
        ccor = res[0]
        pcorr = res[1]
        # if pcorr<0.05:
        #     continue
        if (
            m < minM
        ):  # if small sample size or the correlation between bx and by is ~significant; changed to be 0 bc I actually dont
            continue
        counter = counter + 1
        usedGenes.append(ggs[_])
        ms.append(m)
        ccs.append(ccor)
        ccps.append(pcorr)
        # print(_)
        if counter > nGenesToTry:
            break
        prunedLD, alpha, deco = regularizeLD(
            prunedLD, m + 1, m + 1, 0.1, epsilon=1e-4, toCorr=True, justMCP=True
        )
        Theta = deco[1] @ np.diag(1 / deco[0]) @ deco[1].T
        BTB = np.linalg.inv(bx.T @ Theta @ bx)
        BTa = bx.T @ Theta @ by
        est = BTB @ BTa
        D = np.diag((by - bx @ est).squeeze())
        BT = bx.T @ Theta
        rob = BTB @ BT @ D @ prunedLD @ D @ BT.T @ BTB
        # 6)
        se = np.diag(rob) ** 0.5
        stat = ((est / se) ** 2)[1]  # only care about causal estimate, not intercept
        cstats.append(stat)
    # 7)
    infl = np.median(cstats) / stats.chi2.ppf(0.5, 1) if len(cstats) > 10 else 1
    # finding BP positions of genes used in inflation calculation
    bps = genedf[genedf["Gene"].isin(usedGenes)]["geneBP"].values
    return infl, ccs, ccps, usedGenes, ms, bps
