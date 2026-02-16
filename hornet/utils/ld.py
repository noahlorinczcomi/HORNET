# standard imports
from scipy.stats import norm
import scipy
from numpy import matrix
from scipy import stats
from statsmodels.regression.quantile_regression import QuantReg
import random
import numpy as np
import math

# custom imports
# TODO: fix this import
from MNet import *  # this is a source file in the HORNET directory


def subSamCalcR(Zmat, k, niter=1000, repl=True):
    # Rubins imputation standard errors for rs
    p = Zmat.shape[1]
    n = Zmat.shape[0]
    rs, vBetween, vWithin = np.zeros((p, p)), np.zeros((p, p)), np.zeros((p, p))
    muHat = np.array(Zmat.corr())
    for jj in range(0, niter):
        subSam = Zmat.iloc[:, :].sample(n=k, replace=repl, axis=0)
        cc = np.array(subSam.corr())
        rs = rs + cc
        vWithin = (
            vWithin + ((1 - cc**2) ** 2) / n
        )  # see The Standard Deviation of the Correlation Coefficient, JASA
        vBetween = vBetween + (cc - muHat) ** 2
    rs = rs / niter
    vWithin = vWithin / niter
    vBetween = vBetween / niter
    vTotal = vWithin + vBetween + vBetween / niter
    return (
        rs,
        vTotal,
    )  # correlation matrix estimate, estimate of variance of correlation matrix


# this one subsamples from within evenly-spaced blocks, assuming the data in Zmat are sorted according to BP
def blockSubSamCalcR(ZmatP, k, niter, verbose):
    if verbose == True:
        print("data must be sorted by base pair position")
    p = Zmat.shape[1]
    n = Zmat.shape[0]
    rs, vBetween, vWithin = np.zeros((p, p)), np.zeros((p, p)), np.zeros((p, p))
    blockInds = list(np.linspace(0, Zmat.shape[0] - 1, k))
    blockInds = [int(blockInds[x]) for x in range(0, len(blockInds))]
    muHat = Zmat.loc[blockInds, :].corr()
    for jj in range(0, niter):
        bbInds = [
            random.randint(blockInds[bb], blockInds[bb + 1])
            for bb in range(0, len(blockInds) - 1)
        ]
        subSam = Zmat.iloc[bbInds, :]
        cc = np.array(subSam.corr())
        rs = rs + cc
        vWithin = (
            vWithin + ((1 - cc**2) ** 2) / n
        )  # see The Standard Deviation of the Correlation Coefficient, JASA
        vBetween = vBetween + (cc - muHat) ** 2
    rs = rs / niter
    vWithin = vWithin / niter
    vBetween = vBetween / niter
    vTotal = vWithin + vBetween + vBetween / niter
    return (
        rs,
        vTotal,
    )  # correlation matrix estimate, estimate of variance of correlation matrix


def LDShrinkUpper(ld, tau, snps, missMethod="drop"):
    # if missMethod=='drop':
    #     ld=np.nan_to_num(ld, nan=1)
    lld = ld.copy()
    np.fill_diagonal(lld, 0)
    mmaxs = np.max(lld, axis=0)
    drop = [-1]
    if np.max(mmaxs) < tau:
        return ld, snps
    else:
        if max(mmaxs) > tau:
            drop = [whichMax(mmaxs)]  # finds first occurrence of largest value>tau
        lsnps = list(snps)
        k = 0
        while (-1 in drop) == False:
            k = k + 1
            lld = np.delete(lld, drop[0], axis=0)
            lld = np.delete(lld, drop[0], axis=1)
            # lld=subR(lld,drop);
            del lsnps[drop[0]]
            mmaxs = np.max(lld, axis=0)
            if max(mmaxs) > tau:
                drop = [whichMax(mmaxs)]
            else:
                drop = [-1]
        np.fill_diagonal(lld, 1)
        return lld, lsnps


# MaxDet function for filling in correlation matrix used by MRBEE (requires the function from above: flatten_list())
def MaxDet(R, doToBadSolutions="drop"):
    if (doToBadSolutions == "drop" or doToBadSolutions == "zero") == False:
        raise ValueError("argument doToBadSolutions should be one of drop or zero")
    runningR = R.copy()
    ns = []
    bad_variates = []
    # for each row
    for i_ in range(0, runningR.shape[0]):
        # for each column
        for j_ in range(0, runningR.shape[0]):
            if (
                i_ > j_ or i_ == j_
            ):  # working on upper-triangle only (all will get filled in)
                continue
            else:
                if math.isnan(R[i_, j_]) == False:
                    continue
                else:
                    tempR = R.copy()
                    # reorder, putting target at (n-1)xn and nx(n-1) positions
                    inds_ = list(np.array(range(0, tempR.shape[0])))
                    # (i_, j_) has a missing value
                    # steps
                    # 1) put j_ column last and fix rows
                    # find index location of j_ in inds_ and put j_ at the end
                    a = inds_[: (inds_.index(j_))]
                    b = inds_[(inds_.index(j_) + 1) :]
                    abc = [a, b, j_]
                    abc = flatten_list(abc)
                    tempR = tempR[:, abc][
                        abc, :
                    ]  # does not change the row that missing is on
                    # 2) put i_ row last and fix columns
                    a = inds_[: (inds_.index(i_))]
                    b = inds_[(inds_.index(i_) + 1) :]
                    abc = [a, b, i_]
                    abc = flatten_list(abc)
                    tempR = tempR[abc, :][:, abc]
                    # now missing value is in (n-1)xn and nx(n-1) spot
                    # remove other missing, after making target non-missing
                    n = tempR.shape[0]
                    tempR[n - 1, n - 2] = 99
                    tempR[n - 2, n - 1] = 99
                    # remove any rows/columns with missing
                    # first remove variables where the correlation with the target is missing
                    lastColMiss = np.isnan(tempR[:, -1]) == False
                    tempR = tempR[:, lastColMiss][lastColMiss, :]
                    # now remove other missing
                    # the row with the target should no longer have any missing
                    # make upper triangle 0
                    tempRLower = np.tril(tempR)
                    mask = np.isnan(tempRLower).any(axis=1) == False
                    tempR = tempR[mask, :][:, mask]
                    # impute single missing with 0 and make positive definite (I hope this will ensure solutions stay in {-1,1}
                    # mask=numpy.isnan(tempR)
                    # tempR[mask]=0
                    tempR, alpha, deco = posDefifyCorrMat(
                        tempR, 0.001
                    )  # keep positive definite to avoid bad solutions
                    # identify sub-matrices used for MaxDet solution
                    n = tempR.shape[0]
                    A11 = tempR[: (n - 2), 0 : (n - 2)]
                    B = tempR[(A11.shape[0]), 0 : A11.shape[1]]
                    C = tempR[tempR.shape[0] - 1, 0 : A11.shape[0]]
                    sol = B.T @ np.linalg.inv(A11) @ C
                    ns.append(tempR.shape[0] - 1)
                    if abs(sol) > 1:
                        bad_variates.append((i_, j_))
                        # print(str(i_)+'and'+str(j_)+'sol='+str(sol))
                        if doToBadSolutions == "zero":
                            sol = 0
                    runningR[i_, j_] = sol
                    runningR[j_, i_] = sol
    # It seems like one variable can cause multiple problems and the variable will always
    # be the first value in the tuples which are listed in bad_variates.
    # find the variables whose indices are listed first in each of the tuples in bad_variates
    bad_inds = set([bad_variates[x][0] for x in range(0, len(bad_variates))])
    bad_inds = list(bad_inds)
    goodMask = np.ones((runningR.shape[0],), dtype="bool")
    goodMask[bad_inds] = False
    nDropped = 0 if doToBadSolutions == "zero" else len(bad_inds)
    if doToBadSolutions == "drop":
        runningR = runningR[goodMask, :][:, goodMask]
    return runningR, bad_inds, nDropped


def fullRankMatrixCompletion(
    A, lambda_, missImp=0, thresholding="soft", scad_a=3.7, eps=0.01, max_iter=100
):
    nn = A.shape[0] * A.shape[1]
    nns = []
    k = 0
    Ac = A.copy()
    Omega = np.isnan(Ac)
    Ac[Omega] = (
        missImp  # we hypothesize the mean of missing to be 0 (these are just starting values; the program will add some error)
    )
    lam = lambda_  # may affect results
    while (nn > eps) & (k < max_iter):
        k = k + 1
        deco = np.linalg.svd(Ac, full_matrices=False)
        newd = Sop(deco[1], lam)  # soft-thresholing
        # newd=scad(deco[1],scad_a,lam) # SCAD
        lam = (
            np.min(newd) / 2
        )  # kind of arbitrary, but don't want to set the smallest eigenvalue to 0
        AcPred = deco[0] @ np.diag(newd) @ deco[2]  # no need to transpose the last
        AcPred[Omega == False] = Ac[Omega == False]  # don't change observed values
        nn = nuclearNorm(Ac - AcPred)
        nns.append(nn)
        Ac = AcPred.copy()
    return Ac, k, nns


def taperLD(R, k):
    firstCol = np.linspace(1, R.shape[0], R.shape[0])
    firstCol[:k] = 1
    To = 1 / scipy.linalg.toeplitz(firstCol)
    Rtapered = R * To  # will retain 1's on diagonal
    return Rtapered


def bandLD(R, k):
    if k > R.shape[0]:
        return R
    else:
        firstCol = np.zeros((R.shape[0],))
        firstCol[:k] = 1
        To = scipy.linalg.toeplitz(firstCol)
        Rbanded = R * To
        return Rbanded


def regularizeLD(R, kTaper, kBand, lambda_, epsilon=1e-4, toCorr=True, justMCP=True):
    # epsilon: minimum eigenvalue final matrix is allowed to have
    # toCorr: If True, converts solution to a correlation matrix, else diagonal elements may be outside of {-1,1}
    if justMCP:
        r3 = mcp(R, lambda_)
        r4a, alpha, deco = posDefifyCorrMat(r3, epsilon)
    else:
        r1 = taperLD(R, kTaper)
        r2 = bandLD(r1, kBand)
        r3 = softThresholdLD(r2, lambda_)
        r4a, alpha, deco = posDefifyCorrMat(r3, epsilon)
    if toCorr:
        D = 1 / np.diag(r4a) ** 0.5
        D = np.diag(D)
        r4a = D @ r4a @ D
    return r4a, alpha, deco


def posDefifyCorrMat(ndmatrix, epsilon=1e-4):
    # ndmatrix: pxp square negative-definite matrix of measurement error correlations/covariances
    # since it is a correlation matrix, abs max off-diagonal value is 1
    deco = np.linalg.eig(ndmatrix)
    eig0 = min(np.real(deco[0]))
    eigP = max(np.real(deco[0]))
    if eig0 < epsilon:
        mu = max([epsilon, (eig0 + eigP) / 2])
        alpha = (mu - epsilon) / (
            mu - eig0
        )  # optimal choice given by Choi et al: https://doi.org/10.1016/j.jmva.2018.12.002
        Sigmahat = alpha * ndmatrix + (1 - alpha) * mu * np.eye(
            ndmatrix.shape[0]
        )  # the solution
        deco = np.linalg.eig(Sigmahat)
    else:
        alpha = 1
        Sigmahat = ndmatrix.copy()
    deco = (np.real(deco[0]), np.real(deco[1]))
    return Sigmahat, alpha, deco
