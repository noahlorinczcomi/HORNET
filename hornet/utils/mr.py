# standard imports
from scipy.stats import norm
import warnings
from statsmodels.regression.quantile_regression import QuantReg
import random
import pandas
import numpy as np

# custom imports
from MNet import *  # this is a source file in the HORNET directory


def estMRBEE(Bhat, Ahat, UU, UV, VV, Rinv, subff=True):
    # Bhat,UU: exposure
    # Ahat,VV: outcome
    # Rinv: inverse of LD matrix
    m = len(Bhat[:, 0])
    q = len(Ahat[0])
    p = Bhat.shape[1]
    est = np.linalg.inv(Bhat.T @ Rinv @ Bhat - m * UU) @ (Bhat.T @ Rinv @ Ahat - m * UV)
    F = Bhat.T @ Rinv @ Bhat - m * UU
    res0 = Ahat.squeeze() - (Bhat @ est).squeeze()
    ResMat = np.zeros((len(res0), len(res0)))
    np.fill_diagonal(ResMat, res0.squeeze())
    if subff == True:
        ff = np.zeros((m, p))
        dd = flatten_list([l.tolist() for l in UV - UU @ est])
        for i in range(0, ff.shape[0]):
            ff[i, :] = dd
        E = ResMat @ Bhat - ff
    else:
        E = ResMat @ Bhat
    V = E.T @ Rinv @ E
    Varest = np.linalg.inv(F) @ V @ np.linalg.inv(F)
    return est, Varest


def bootMRBEE(Bhat, Ahat, UU, UV, VV, R, niter=500):
    # R: LD matrix
    m = len(Bhat[:, 0])
    p = Bhat.shape[1]
    ests = np.zeros((niter, p))
    for j in range(0, niter):
        ind = random.choices(range(0, m), k=m)
        kBhat = Bhat[ind, :]
        kAhat = Ahat[ind, :]
        kR1 = R[ind, :]
        kR = kR1[:, ind]
        kR = np.linalg.pinv(
            kR
        )  # generalized inverse because not guaranteed to be posdef
        est, V = estMRBEE(kBhat, kAhat, UU, UV, VV, kR)
        estl = [l.tolist() for l in est]
        ests[j, :] = flatten_list(estl)
    cc = pd.DataFrame(ests).cov()
    ests = np.array([np.mean(ests[:, x]) for x in range(0, bx.shape[1])])
    cc = np.array(cc)
    return ests, cc


def imrbee(
    Bhat,
    Ahat,
    UU,
    UV,
    VV,
    Rinv,
    R,
    PleioPThreshold,
    boot=False,
    niter=1000,
    initial="bee",
    max_iter=15,
):
    warnings.filterwarnings("ignore")
    # R: LD matrix
    # Rinv: inverse of LD matrix
    m = len(Bhat[:, 0])
    p = len(Bhat[0, :]) + 1
    q = len(Ahat[0, :])  # +1 on p bc I will add an intercept without asking
    Bhat_ = np.column_stack(([1] * Bhat.shape[0], Bhat))
    UU_ = np.column_stack(([0] * (p - 1), UU))
    UU_ = np.row_stack(([0] * p, UU_))
    UV_ = np.row_stack((0, UV))
    if initial == "robust":
        mod = QuantReg(Ahat, Bhat_)  # includes an intercept
        res = mod.fit(q=0.5)
        est0 = res.params.reshape((p, 1))
        V0 = res.cov_params()
        pleioPs0 = SpleioP(Bhat_, Ahat, UU_, UV_, VV, est0, V0)
    else:
        est0, V0 = estMRBEE(Bhat_, Ahat, UU_, UV_, VV, Rinv)
        pleioPs0 = SpleioP(Bhat_, Ahat, UU_, UV_, VV, est0, V0)
    k = 0
    thetadiff = 1
    diffs = []
    kR = Rinv.copy()
    mask = np.ones((m,), dtype="bool")
    Outliers0 = [
        i for i in range(0, len(pleioPs0)) if pleioPs0[i] < 0.05
    ]  # the initial threshold can be stricter than the one used in iteration
    mask[Outliers0] = False
    Outliers0 = np.array(Outliers0)
    Outliers = Outliers0.copy().tolist()
    while (
        (k < max_iter)
        & (thetadiff > (0.0001 * (p * q + p)))
        & (len(Outliers) < (m - 10))
    ):
        k = k + 1
        kBhat = Bhat_[mask]
        kAhat = Ahat[mask]
        Rsub = R[mask, :][:, mask]
        kR = np.linalg.pinv(Rsub)
        fitkEst, fitkVar = estMRBEE(kBhat, kAhat, UU_, UV_, VV, kR)
        pleioPs0 = SpleioP(Bhat_, Ahat, UU_, UV_, VV, fitkEst, fitkVar)
        Outliers = [i for i in range(0, len(pleioPs0)) if pleioPs0[i] < PleioPThreshold]
        mask[Outliers] = 0
        if len(Outliers) == kBhat.shape[0]:
            thetadiff = 0
            est0 = np.array([1] * p) * float("nan")
            V0 = np.zeros((p, p)) * float("nan")
            Outliers = list(range(0, m))
            kR = float("nan")
        elif len(Outliers) == 0:
            thetadiff = 0
            est0, V0, Outliers, kR = fitkEst, fitkVar, Outliers0, kR
        else:
            thetadiff = np.sum(abs(fitkEst - est0))
            est0 = fitkEst
            V0 = fitkVar
            diffs.append(thetadiff)
    if (boot is True) & (len(Outliers) < m):
        if len(Outliers) == 0:
            kBhat = Bhat.copy()
            kAhat = Ahat.copy()
            Rsub = R.copy()
            kR = np.linalg.inv(Rsub)
        else:
            kBhat = np.delete(Bhat, (Outliers), axis=0)
            kAhat = np.delete(Ahat, (Outliers), axis=0)
            Rsub = subR(R, Outliers)
            kR = np.linalg.pinv(Rsub)
        est0, V0 = bootMRBEE(kBhat, kAhat, UU, UV, VV, kR, niter=niter)
    warnings.filterwarnings("default")
    return est0, V0, Outliers, kR, k


def IVW(bx, by, Theta, ogZs, uniIVPThreshold=5e-8, has_intercept=False):
    # ogZs: original Z-stats (in the case of a transformation of bx by LD structure)
    # Theta: inverse of LD matrix
    # this will perform univariable and multivariable MR
    # uniIVPThreshold: if the user wants to restrict univariable MR IVs to only those significant at some threshold (recall bx and by are Z-scores)
    m = bx.shape[0]
    p = bx.shape[1]
    if has_intercept == False:
        bx_ = np.column_stack(([1] * m, bx))
        # MVMR
    BTB = bx_.T @ Theta @ bx_
    hasInverse = np.linalg.matrix_rank(BTB) == BTB.shape[0]
    if hasInverse:
        BTBinv = np.linalg.inv(BTB)
    else:
        BTBinv = np.linalg.pinv(BTB)
    mvmrest = BTBinv @ bx_.T @ Theta @ by
    mvmrest = mvmrest.reshape((p + 1,))
    res = by.reshape((m,)) - bx_ @ mvmrest
    sig2 = sum(res**2) / (m - p - 1)
    mvmrse = np.diag(sig2 * BTBinv) ** 0.5
    mvmrest = mvmrest[1:]  # take off intercept-relevant terms
    mvmrse = mvmrse[1:]  # take off intercept-relevant terms
    # univariable
    UniEst = []
    UniSE = []
    q0 = norm.ppf(1 - uniIVPThreshold)
    for _ in range(0, p):
        bxj_ = bx[:, _]
        bxj_ = np.column_stack(([1] * m, bxj_))
        mask = np.abs(ogZs[:, _]) > q0
        newm = sum(mask)
        if (
            newm < 5
        ):  # if <5 SNPs meet P-value threshold (some genes may have SNPs like this since we keep all SNPs in MRJones with P<tau in a chi-square joint test)
            UniEst.append(float("nan"))
            UniSE.append(float("nan"))
        else:
            bxj_ = bxj_[mask]
            by_ = by[mask]
            Theta_ = Theta[mask, :][:, mask]
            btbj = np.linalg.inv(bxj_.T @ Theta_ @ bxj_)
            uniest = btbj @ bxj_.T @ Theta_ @ by_
            uniest = uniest.reshape((2,))
            res = by_.reshape((newm,)) - bxj_ @ uniest
            sig2 = sum(res**2) / (newm - 2)
            unise = np.diag(sig2 * btbj) ** 0.5
            UniEst.append(uniest[1])
            UniSE.append(unise[1])
    dd = {
        "IVW_MVMR_Est": mvmrest,
        "IVW_MVMR_SE": mvmrse,
        "IVW_UNIMR_Est": UniEst,
        "IVW_UNIMR_SE": UniSE,
    }
    pp = pd.DataFrame(dd)
    return pp


def MRBEERidge(Bhat, Ahat, UU, UV, VV, Rinv, nLambda=100):
    # estimator is simple; lambdaMax determined automatically
    m = Bhat.shape[0]
    p = Bhat.shape[1]
    lambdaMax = np.max(abs(Bhat.T @ Rinv @ Ahat) / m) ** (1 / 2)
    lls = np.linspace(0, lambdaMax, nLambda)
    logmse = []
    for i in range(0, len(lls)):
        ee = np.linalg.inv(Bhat.T @ Rinv @ Bhat - m * UU + lls[i] * np.eye(p)) @ (
            Bhat.T @ Rinv @ Ahat - m * UV
        )
        res = Ahat - Bhat @ ee
        toAdd = np.log(sum(res**2)[0])
        logmse.append(toAdd)
        if i > 1:
            if logmse[-1] > logmse[-2]:
                break
    # final estimation
    lambda0 = lls[i - 1]
    Ip = makeI(p)
    ee0 = np.linalg.inv(Bhat.T @ Rinv @ Bhat - m * UU + lambda0 * Ip) @ (
        Bhat.T @ Rinv @ Ahat - m * UV
    )
    # variance-covariance
    Fhat = (1 / m) * (Bhat.T @ Rinv @ Bhat + lambda0 * Ip - m * UU)
    ResMat = np.zeros((m, m))
    np.fill_diagonal(ResMat, flatten_list(Ahat - Bhat @ ee0))
    ff = np.zeros((m, p))
    dd = flatten_list([l.tolist() for l in UV - UU @ ee0 + lambda0 * ee0])
    for i in range(0, ff.shape[0]):
        ff[i, :] = dd
    E = ResMat @ Bhat - ff
    Vhat = E.T @ Rinv @ E / m
    Varhat = 1 / m * np.linalg.inv(Fhat) @ Vhat @ np.linalg.inv(Fhat)
    return ee0, Varhat, lambda0


def MRJones(
    by,
    bX,
    LD,
    Rxx,
    rxy,
    lamvec=np.linspace(0.05, 0.3, 10),
    tauvec=np.linspace(5, 25, 10),
    rho_theta=1,
    rho_gamma=1,
    max_iter=30,
    max_eps=0.005,
    intercept=True,
    WIB=True,
    WIB_thres=5,
    bicfactor=1,
    normmax=3,
):
    warnings.filterwarnings("ignore")
    n = len(by.squeeze())
    if intercept:
        bX = np.column_stack(([1] * n, bX))
        Rxx = np.column_stack(([0] * Rxx.shape[0], Rxx))
        Rxx = np.row_stack(([0] * (Rxx.shape[1]), Rxx))
        rxy = np.row_stack(([0], rxy))
    p = bX.shape[1]
    q = len(lamvec)
    w = len(tauvec)
    Theta = np.linalg.inv(LD)
    TC = np.linalg.cholesky(LD).T
    BT = bX.T @ Theta
    BTB = BT @ bX
    BTa = BT @ by
    adjvec = np.diag(BTB - n * Rxx)
    adjvec.flags.writeable = True
    mask = adjvec < 0
    adjvec[mask] = 1 / (n**3)
    adjvec = adjvec**0.5
    adjvec = adjvec / np.max(adjvec)

    if WIB:
        ss = np.diag(BTB) / n
        mask = ss < WIB_thres
        qq = np.diag(Rxx)
        qq.flags.writeable = True
        qq[mask] = (qq[mask] * ss[mask] / WIB_thres) ** 0.5
        Rxx = np.diag(qq) @ Rxx @ np.diag(qq)

    Thetarho = np.linalg.inv(LD + rho_gamma * np.diag([1] * n))
    by1 = (TC @ by).squeeze()
    if intercept:
        bX1 = TC @ bX[:, 1:]
        fit0 = QuantReg(by, bX1)
        res = fit0.fit(
            q=0.5
        )  # can give a warning if max iterations reached... how to silence?
        theta_ini = res.params
        theta_ini = theta_ini.reshape((theta_ini.shape[0], 1))
        theta_ini.flags.writeable = True
        mask = np.abs(theta_ini) < 0.001
        theta_ini[mask] = 0
        theta_ini = np.row_stack((0, theta_ini))
        adjvec[0] = 1
    else:
        bX1 = TC @ bX
        fit0 = QuantReg(by, bX1)
        res = fit0.fit(q=0.5)
        theta_ini = res.params
        theta_ini = theta_ini.reshape((theta_ini.shape[0], 1))
        theta_ini.flags.writeable = True
        mask = np.abs(theta_ini) < 0.001
        theta_ini[mask] = 0
        theta_ini = np.row_stack((0, theta_ini))

    normmax = max(np.sqrt(sum(theta_ini**2)) * normmax)
    normmax = normmax if normmax > 1 else 1
    r = (by - bX @ theta_ini).squeeze()
    gamma_ini = (Theta @ r).squeeze()
    mask = abs(gamma_ini / 3) < mad(gamma_ini)
    gamma_ini[mask] = 0

    lamvec = np.sort(lamvec)
    tauvec = np.sort(tauvec)
    Btheta = np.zeros((p, q, w))
    Bgamma = np.zeros((n, q, w))
    Bbic = np.zeros((q, w))

    for ss in range(0, q):
        for sss in reversed(range(0, w)):
            theta = theta_ini.copy().squeeze()
            theta1 = theta.copy().squeeze()
            gamma = gamma_ini.copy().squeeze()
            gamma1 = gamma.copy().squeeze()
            u1 = rho_theta * (theta - theta1)
            u1 = u1.squeeze()
            u2 = rho_gamma * (gamma - gamma1)
            u2 = u2.squeeze()
            theta2 = theta1 * 0
            error = 1
            iter_ = 0
            while (error > max_eps) & (iter_ < max_iter):
                iter_ = iter_ + 1
                theta2 = theta1.copy()
                effn = sum(gamma1 == 0)
                Hinv = np.linalg.inv(BTB - effn * Rxx + rho_theta * np.diag([1] * p))
                theta = Hinv @ (
                    BTa.squeeze()
                    - bX.T @ gamma
                    - (rxy.squeeze()) * effn
                    - (u1.squeeze() - rho_theta * theta1.squeeze())
                )
                if np.sqrt(sum(theta.squeeze() ** 2)) > normmax:
                    theta = theta / np.sqrt(sum(theta**2)) * normmax

                theta1 = mcp(
                    theta + u1.squeeze() / rho_theta,
                    lamvec[ss] / adjvec / rho_theta,
                    ga=3,
                )
                gamma = Thetarho @ (
                    by.squeeze() - bX @ theta - u2.squeeze() + rho_gamma * gamma1
                )
                gamma1 = mcp(
                    gamma + u2.squeeze() / rho_gamma, tauvec[sss] / rho_gamma, ga=3
                )
                u1 = u1 + rho_theta * (theta - theta1)
                u2 = u2 + rho_gamma * (gamma - gamma1)
                error = (
                    sum((theta2.squeeze() - theta1.squeeze()) ** 2) if iter_ > 3 else 1
                )
                # end while

            Btheta[:, ss, sss] = theta1.squeeze()
            Bgamma[:, ss, sss] = gamma1.squeeze()
            df1 = sum(theta1 != 0)
            df2 = sum(gamma1 != 0)
            r = by.squeeze() - bX @ theta1 - LD @ gamma1
            varr = mad(TC @ r) ** 2
            if varr == 0:  # if MRJones wanted all to be pleiotropy
                rss = 1e5  # something very large bc I never want this solution
            else:
                rss = sum(r * (Theta @ r) / varr)
            Bbic[ss, sss] = (
                n * np.log(rss)
                + (np.log(n) + np.log(p) * bicfactor) * df1
                + (np.log(n) + np.log(n) * bicfactor) * df2
            )
            # end for sss
    # end for ss
    istar, jstar = np.where(Bbic == Bbic.min())
    istar = istar[0]
    jstar = jstar[0]  # may be multiple minimums? Choose first
    theta = Btheta[:, istar, jstar]
    gamma = Bgamma[:, istar, jstar]
    # ignore everything intercept-related from the beginning
    if intercept:
        theta[0] = 0  # instead of deleting bc I use bX[:,masktheta] below
        covg = np.zeros((p - 1, p - 1))
    else:
        covg = np.zeros((p, p))
    indtheta = [i for i in range(0, len(theta)) if theta[i] != 0]
    indgamma = [i for i in range(0, len(gamma)) if gamma[i] != 0]
    effn = n - len(indgamma)
    masktheta = np.zeros((len(theta)), dtype="bool")
    masktheta[indtheta] = True
    maskgamma = np.ones((len(r)), dtype="bool")
    maskgamma[indgamma] = False

    if (sum(indtheta) > 0) & (
        sum(maskgamma) > 2
    ):  # must have at least 2 non pleiotropic IVs
        Z = bX[:, masktheta][maskgamma, :]
        R0inv = np.linalg.inv(LD[maskgamma, :][:, maskgamma])
        r = by.squeeze() - bX @ theta - LD @ gamma.squeeze()
        r = r[maskgamma]
        D = np.diag(r.squeeze())
        V = Z.T @ R0inv @ D @ (LD[maskgamma, :][:, maskgamma]) @ D @ R0inv @ Z
        Hinv = np.linalg.inv(
            Z.T @ R0inv @ Z - sum(maskgamma) * Rxx[masktheta, :][:, masktheta]
        )
        covg = Hinv @ V @ Hinv

    # if intercept:
    #     intercept=theta[0]
    #     theta=theta[1:]
    #     # covg=covg[1:,:][:,1:] # don't need this bc I already made theta for intercept 0 above
    #     bX=bX[:,1:]
    if intercept:
        theta = theta[1:]

    dout = {}
    dout["theta"] = theta
    dout["gamma"] = gamma
    dout["selectedgenes"] = indtheta
    dout["pleiotropyIVs"] = indgamma
    dout["covg"] = covg
    dout["Bic"] = Bbic
    dout["Btheta"] = Btheta
    dout["Bgamma"] = Bgamma
    dout["intercept"] = intercept
    warnings.filterwarnings("default")
    return dout


def singleSNP(bx, by, cn):  # assumes Z-stat standardization
    m = bx.shape[0]
    p = bx.shape[1]
    by = by.squeeze()
    wald = []
    waldse = []
    for _ in range(0, p):
        bxp = bx[:, _].squeeze()
        maxix = np.argmax(abs(bxp))
        wald.append(by[maxix] / bxp[maxix])
        waldse.append(1 / bxp[maxix] ** 2)
    wald = np.array(wald)
    waldse = np.array(waldse)
    ps = [2 * norm.cdf(-abs(wald[_] / waldse[_]), 0, 1) for _ in range(0, p)]
    dfout = pandas.DataFrame({"Gene": cn, "SingleSNPEst": wald, "SingleSNP_P": ps})
    return dfout


def MRBEEHuber(
    bx, by, R, UU, UV, gamma=0.5, initial="ivw", eps=0.001, max_iter=15, boot=False
):
    m = bx.shape[0]
    p = bx.shape[1]
    # take away intercept
    by = by - np.mean(by)
    for _ in range(0, bx.shape[1]):
        bx[:, _] = bx[:, _] - np.mean(bx[:, _])
    w, v = np.linalg.eig(R)
    Rinvsqrt = v @ np.diag(1 / w**0.5) @ v.T
    bx_ = Rinvsqrt @ bx  # transforming up front
    by_ = Rinvsqrt @ by  # transforming up front
    if np.all(bx[:, 0] == 0):  # if user gave bx with an intercept
        bx_[:0] = 1  # ensure it stays like it should
    if initial == "robust":
        mod = QuantReg(by_, bx_)  # includes an intercept
        res = mod.fit(q=0.5)
        est0 = res.params.reshape((p, 1))
        res0 = by_ - bx_ @ est0
    else:
        est0 = np.linalg.inv(bx_.T @ bx_) @ bx_.T @ by_
        res0 = by_ - bx_ @ est0
    error = 1 + eps
    k = 0
    while (error > eps) & (k < max_iter):
        k = k + 1
        w0 = huberWeight(res0, gamma)
        D = np.diag(w0.squeeze())
        # w0=w0/numpy.max(w0)
        adjfactor = sum(w0)
        H = bx_.T @ D @ bx_ - adjfactor * UU
        kk = 0.95
        while (ispd(H) == False) & (kk > 0):
            kk = kk - 0.05
            UU = UU * kk + (1 - kk) * np.eye(p)
        estk = np.linalg.inv(bx_.T @ D @ bx_ - adjfactor * UU) @ bx_.T @ D @ by_
        error = (est0 - estk).T @ (est0 - estk)
        est0 = estk.copy()
        res0 = by_ - bx_ @ est0
    adjfactor = sum(w0)
    if boot:
        varest = bootVarWLS(bx_, by_, w0, UU)
    else:
        s2 = sum((w0 * res0) ** 2) / (
            adjfactor - 2
        )  # this may underestimate variance (assumes fixed weights)
        Hinv = np.linalg.inv(bx_.T @ D @ bx_ - adjfactor * UU)
        XWX = bx_.T @ (D**2) @ bx_
        varest = s2 * Hinv @ XWX.T @ Hinv
    return est0, varest, w0
