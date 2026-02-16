# standard imports
from scipy.stats import norm
from scipy import stats
from statsmodels.regression.quantile_regression import QuantReg
import pandas
import numpy as np

# custom imports
from MNet import *  # this is a source file in the HORNET directory


def perGeneSMRBias(
    bx,
    by,
    ld0,
    sparsePrecision,
    ogZs,
    UU,
    UV,
    cn,
    IVPThresh=5e-5,
    r2Thresh=0.25,
    minNIVs=10,
):
    # ogZs: original Z-stats in case a transformation of bx was made
    # this function uses gene groups defined by the MVMR code. For each gene, it estimates its total causal effect unconditional on
    # the other genes in the group. It then estimates the extent to which this causal estimate is biased by CHP.
    # this therefore doesn't require any additional LD calculations
    # mrbee estimates
    p = bx.shape[1]
    m = bx.shape[0]
    VV = np.array([1]).reshape((1, 1))
    e0, v0, o0, k0, kiter0 = imrbee(
        bx, by, UU, UV, VV, sparsePrecision, ld0, PleioPThreshold=0.05 / m, boot=False
    )  # strict-ish PleioPThreshold
    e0 = e0[1:]
    v0 = v0[1:, :][:, 1:]  # take off intercept-relevant terms
    mask0 = numpy.ones((bx.shape[0],), dtype="bool")
    if len(o0) > 0:
        mask0[o0] = False
    _0 = ld0[mask0, :][:, mask0]
    _0sq = np.linalg.cholesky(_0).T
    _0 = np.linalg.pinv(_0)
    A1 = (
        numpy.linalg.pinv(bx[mask0, :].T @ _0 @ bx[mask0, :] - m * UU)
        @ bx[mask0, :].T
        @ _0
    )
    res0 = by.reshape((by.shape[0],)) - bx @ e0.reshape((e0.shape[0],))
    res0 = res0[mask0]
    mjs = []
    jointStats = []
    jointPs = []
    pleioMeans = []
    pleioMeansP = []
    fronorm = []
    diffStats = []
    diffStatPs = []
    uniests = []
    uniestses = []
    q0 = norm.ppf(1 - IVPThresh)
    for _ in range(0, p):
        inds0 = np.array(range(0, bx.shape[0]))
        bxj = bx[:, _].reshape((m, 1))
        mask = np.abs(ogZs[:, _].reshape((m, 1))) > q0
        mask = mask.reshape((m,))
        R11 = ld0[mask, :][:, mask]
        R12 = ld0[mask, :][:, mask == False]
        bxj = bxj[mask, :]
        mj = bxj.shape[0]
        byj = by[mask, :]
        inds0 = inds0[mask]
        if (
            R12.shape[1] == 0
        ):  # if all SNPs in bx are IVs for this gene (ie no other SNPs to compare these ones to)
            mjs.append(float("nan"))
            jointStats.append(float("nan"))
            jointPs.append(float("nan"))
            pleioMeans.append(float("nan"))
            pleioMeansP.append(float("nan"))
            diffStats.append(float("nan"))
            diffStatPs.append(float("nan"))
            uniests.append(float("nan"))
            uniestses.append(float("nan"))
            continue
        # to make the MRBEE total causal effect estimate robust, I will remove SNPs that are in high LD with other SNPs in this region
        # recall that R12 has rows as this gene's IVs and columns as all OTHER gene's IVs
        mask2 = [
            np.max((R12[x, :]) ** 2) < r2Thresh for x in range(0, R12.shape[0])
        ]  # maks where IV not in LD with non-gene SNPs beyond r2>r2Thresh
        mask2 = numpy.array(mask2)
        # if not enough data, cannot estimate anything reliably
        if sum(mask2) < minNIVs:
            if R12.shape[0] == 0:
                fronorm.append(float("nan"))
            else:
                R11inv = numpy.linalg.pinv(R11)
                fronorm.append(
                    numpy.linalg.norm(R12.T @ R11inv @ R11inv @ R12, ord="fro")
                )
            mjs.append(float("nan"))
            jointStats.append(float("nan"))
            jointPs.append(float("nan"))
            pleioMeans.append(float("nan"))
            pleioMeansP.append(float("nan"))
            diffStats.append(float("nan"))
            diffStatPs.append(float("nan"))
            uniests.append(float("nan"))
            uniestses.append(float("nan"))
            continue
        bxj = bxj[mask2, :]
        mj = bxj.shape[0]
        byj = byj[mask2, :]
        UUj = numpy.array([1]).reshape((1, 1))
        UVj = UV[_, :].reshape((1, 1))
        R11 = R11[mask2, :][:, mask2]
        R12 = R12[mask2, :]
        inds0 = inds0[mask2]
        # ld0, and thus R11, will always be positive definite because I already made it so in the code before this function is executed
        R11inv = np.linalg.pinv(R11)
        precisionjsqrt = np.linalg.cholesky(R11inv).T
        mjs.append(bxj.shape[0])
        # need to estimate total causal effect without bias (cannot use MVMR estimates - those are not total)
        # bootsel=True if bxj.shape[0]<30 else False
        bootsel = False  # do not allow bootstrapping - only causes problems
        pleiot = 0.05 / bxj.shape[0]
        est0, V0, Outliers, kR, k = imrbee(
            bxj,
            byj,
            UUj,
            UVj,
            VV,
            R11inv,
            R11,
            PleioPThreshold=pleiot,
            boot=bootsel,
            niter=500,
        )  # conservative PleioPThreshold
        est0 = est0[1:]
        V0 = V0[1:, :][:, 1:]  # take off intercept-relevant terms
        uniests.append(float(est0))
        uniestses.append(float(V0) ** 0.5)
        # Outliers are identified
        pleioPs = SpleioP(bxj, byj, UUj, UVj, VV, est0, V0)
        pleio = [(byj[x] - est0[0] * bxj[x, 0])[0] for x in range(0, mj)]
        pleioPs = np.array(pleioPs)
        pleio = np.array(pleio)
        newmask = np.ones((bxj.shape[0],), dtype="bool")
        # potential problem here ... if all are outliers in IMRBEE (i.e., Outliers=list(range(0,m)), then I can't continue
        if len(Outliers) > 0:
            newmask[Outliers] = False
        if sum(newmask) == 0:  # if all outliers
            quickAppend(
                [jointStats, jointPs, pleioMeans, pleioMeansP, diffStats, diffStatPs]
            )
            continue
        inds0 = inds0[newmask]
        # joint test for no pleio
        V = (
            R11
            + float(est0**2) * R11
            + float(V0) * bxj @ bxj.T
            - 2 * float(est0 * UVj) * R11
        )
        Vinv = numpy.linalg.pinv(V)
        pleio = pleio.reshape((pleio.shape[0], 1))
        chistat = numpy.sum(numpy.diag(Vinv @ pleio @ pleio.T))
        jointStats.append(chistat)
        jointPs.append(1 - stats.chi2.cdf(chistat, pleio.shape[0]))
        # test for unbalanced HP
        jm = numpy.ones((bxj.shape[0],))
        Xi = R11 * (numpy.eye(R11.shape[0]) + float(est0 * UUj) - 2 * float(est0 * UVj))
        eta = 1 / len(pleio) ** 2 * (jm.T @ Xi @ jm)
        stat = 1 / bxj.shape[0] * (jm.T @ precisionjsqrt @ pleio)
        stat = float(stat)
        statP = 2 * norm.cdf(-abs(stat), 0, 1 / bxj.shape[0] ** 0.5)
        pleioMeans.append(float(stat))
        pleioMeansP.append(float(statP))
        # test for change in uni MR est vs MVMR est
        res1 = pleio.copy().reshape((pleio.shape[0],))
        res1 = res1[newmask]  # see above (non-pleio SNPs)
        # need Cov(y1|x1,y2|x2) after removing pleio
        _2 = R11[newmask, :][:, newmask]
        _2sq = np.linalg.cholesky(_2).T
        _2 = np.linalg.pinv(_2)
        A2 = (
            numpy.linalg.pinv(
                bxj[newmask, :].T @ _2 @ bxj[newmask, :] - sum(newmask) * UUj
            )
            @ bxj[newmask, :].T
            @ _2
        )
        R01 = ld0[mask0, :][
            :, inds0
        ]  # mask0: inds of nonpleio SNPs from MVMRBEE; inds0: inds of nonpleio SNPs from UniMRBEE
        Cov_ = numpy.diag(res0) @ R01 @ numpy.diag(res1)
        Cov_ = A1 @ Cov_ @ A2.T
        Cov_ = Cov_[
            _
        ]  # only selecting the covariance corresponding to this gene/exposure (don't care about Cov(thisGeneThetahat, otherGeneThetaHat)
        if v0.shape[0] == 1:
            diff_ = e0[_] - est0
            estvar = V0 + v0[_] - 2 * Cov_
        else:
            diff_ = e0[_] - est0
            estvar = V0 + v0[_, _] - 2 * Cov_
        if estvar < 0:  # don't know why this happens but sometimes seems possible
            diffStats.append(float("nan"))
            diffStatPs.append(float("nan"))
            continue
        diffStat = diff_ / (estvar**0.5)
        diffStatP = 2 * norm.cdf(-abs(diffStat), 0, 1)
        diffStats.append(float(diffStat))
        diffStatPs.append(float(diffStatP))
    data = {
        "gene": cn,
        "UniMRBEEEsts": uniests,
        "UniMRBEESEs": uniestses,
        "MRBEEEsts": flatten_list(e0.tolist()),
        "MRBEESEs": (numpy.diag(v0) ** 0.5).tolist(),
        "SMRm": mjs,
        "SMRpleioJointStats": jointStats,
        "SMRpleioJointPs": jointPs,
        "SMRHPmeans": pleioMeans,
        "SMRHPmeansP": pleioMeansP,
        "OVBiasTestStat": diffStats,
        "OVBiasTestP": diffStatPs,
    }
    df = pandas.DataFrame(data)
    # some notes on the output
    # 'NaN' in all but last column where number of IVs<5; 'NaN' in last col if number of IVs was 0
    return df


#  function to calculate conditional F-statistics
# !!! at the moment, this function may be extremely memory-intensive (consider the kronecker products)
def conditionalFMRBEE(
    bx,
    UU,
    sparsePrecisionSqrt,
    tauLowerQuant=0,
    tauUpperQuant=0,
    opAlphaVal=0,
    method="mrbee",
):
    # MRJones and MRBEE with IMRP adjustment give slightly different estimates of the conditional Fs
    # I trust MRJones more because it also controls for high correlations between exposures, which MRBEE cannot explicitly do
    bx_ = sparsePrecisionSqrt @ bx.copy()
    Theta = makeI(bx.shape[0])
    Fs = []
    for _ in range(0, bx.shape[1]):
        x_ = numpy.delete(bx_, (_), axis=1)
        x_ = x_.reshape((x_.shape[0], bx.shape[1] - 1))
        y_ = bx_[:, _]
        y_ = y_.reshape((y_.shape[0], 1))
        UU_ = numpy.delete(UU, (_), axis=0)
        UU_ = numpy.delete(UU_, (_), axis=1)
        UV_ = UU[:, _]
        UV_ = numpy.delete(UV_, (_), axis=0)
        UV_ = UV_.reshape((UV_.shape[0], 1))
        # if user just wants to use MRBEE (basically) to find conditional F-statistics (much faster)
        if method.lower() == "mrbee":
            ### using MRBEE-IMRP
            fit = imrbee(
                x_,
                y_,
                UU_,
                UV_,
                numpy.array([1]).reshape((1, 1)),
                Theta,
                numpy.linalg.inv(Theta),
                PleioPThreshold=0.05 / (x_.shape[0] ** 0.5),
                boot=False,
            )
            finalEsts = fit[0]  # est0, V0, Outliers, kR, k
            finalEsts = finalEsts[1:]  # ignore intercept
            mask = numpy.ones((finalEsts.shape[0],), dtype="bool")
            estsVars = fit[1][1:, :][
                :, 1:
            ]  # remove items corresponding to intercept (automatically added by MRBEE)
            df1 = sum(mask)
            df2 = 0
        else:
            lamMax = (numpy.max(abs(x_.T @ y_)) / bx.shape[0]) ** 0.5
            lamvec = numpy.linspace(0, lamMax, 20)
            e0 = numpy.linalg.inv(x_.T @ x_ - x_.shape[0] * UU_) @ (
                x_.T @ y_ - x_.shape[0] * UV_
            )  # initial unbiased estimate with no shrinkage or pleiotropy adjustment
            res = y_.reshape((y_.shape[0],)) - (x_ @ e0).reshape((y_.shape[0],))
            kt1 = max(abs(res)) * tauLowerQuant
            kt2 = max(abs(res)) * tauUpperQuant
            tauvec = numpy.linspace(kt1, kt2, 20)
            ### using MRJones
            fit = MRJones(
                x_,
                y_,
                Theta,
                UU_,
                UV_,
                lamvec,
                tauvec,
                rho_theta=1,
                rho_delta=1,
                max_iter=10,
                max_eps=0.01,
            )
            outSearched = BICgridSearch(
                fit, lamvec, tauvec
            )  # finding optimal hyperparameters
            finalEsts = fit["thetamatrix"][
                :, outSearched["i"], outSearched["j"]
            ]  # extracting final causal estimates from output
            mask = finalEsts != 0
            finalDeltas = fit["deltamatrix"][
                :, outSearched["i"], outSearched["j"]
            ]  # extracting final xi/horizontal pleiotropy estiamtes from output
            finalLambda = outSearched["minLambda1"]
            finalTau = outSearched[
                "minLambda2"
            ]  # definining optimal lambda1 and lambda2/tau
            estsVars = MRJonesVariance(
                x_,
                y_,
                finalDeltas,
                Theta,
                UU_,
                UV_,
                finalEsts,
                finalDeltas != 0,
                finalLambda,
                finalTau,
                opAlphaVal,
            )  # estimating variance of causal estimates
            df1 = sum(mask)
            df2 = sum(finalDeltas != 0)

        # chi-square statistic (since I pre-multiplied by sparsePrecisionSqrt, all IVs are independent)
        finalEsts = finalEsts.reshape((finalEsts.shape[0], 1))
        var_ = (
            1
            + float(finalEsts.T @ finalEsts)
            + numpy.diag(x_[:, mask] @ estsVars @ x_[:, mask].T)
        )
        res = y_.squeeze() - (x_[:, mask] @ finalEsts[mask]).squeeze()
        chi = sum(res**2 / var_)
        condF = chi / (x_.shape[0] - df1 - df2 - 1)
        Fs.append(condF)
    return Fs


def SpleioP(Bhat, Ahat, UU, UV, VV, Thetahat, VThetahat):
    m = len(Bhat[:, 0])
    q = len(Ahat[0, :])
    I = makeI(q)
    Spleio = []
    for i in range(0, m):
        K = np.kron(I, Bhat[i, :])
        num = Ahat[i, :] - Thetahat.T @ Bhat[i, :]
        den = (
            VV + Thetahat.T @ UU @ Thetahat + K @ VThetahat @ K.T - 2 * Thetahat.T @ UV
        )
        S = num @ np.linalg.inv(den) @ num.T
        S = float([S][0])
        Spleio.append(S)
    P = [1 - stats.chi2.cdf(Spleio[x], q) for x in range(0, m)]
    return P
