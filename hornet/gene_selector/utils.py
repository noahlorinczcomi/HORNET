# standard imports
import numpy as np


# Gene selector (my alternative to MR-Jones)
def admmHuber(
    bx, by, SigmaUU, gamma, lam, rho=1, eps=1e-2, max_iter=15, stdize_weights=False
):
    # this function expects bx and by to already be transformed out of LD
    m = bx.shape[0]
    p = bx.shape[1]
    by = by.squeeze()
    UV = np.zeros((p, 1))
    thetaT, v0, w0 = MRBEEHuber(
        bx,
        by,
        np.eye(m),
        SigmaUU,
        UV,
        gamma=gamma,
        initial="bee",
        eps=eps,
        max_iter=max_iter,
        boot=False,
    )
    if stdize_weights:
        w0 = w0 / np.max(w0)

    D = np.diag(w0.squeeze())
    adj = sum(w0.squeeze())
    I = np.eye(p)
    epsie = [1]
    k = 0
    muT = np.zeros((p,))  # fixed start
    ts = np.zeros((max_iter, p))
    ws = np.zeros((max_iter, m))
    while (epsie[-1] > eps) & (k < max_iter):
        k = k + 1
        BTB = bx.T @ D @ bx + rho * I - adj * SigmaUU
        BTa = bx.T @ D @ by - muT
        betaTp1 = np.linalg.inv(BTB) @ BTa
        betaTp1 = betaTp1.reshape(
            p,
        )
        thetaTp1 = scad(betaTp1 + muT / rho, lam / rho, 3.7)
        muTp1 = muT + rho * (betaTp1 - thetaTp1)
        oops = np.sum((betaTp1 - thetaTp1) ** 2)
        epsie.append(oops)
        thetaT = thetaTp1.reshape(
            p,
        )
        muT = muTp1
        res = by.squeeze() - bx @ (
            thetaT.reshape(
                p,
            )
        )
        w0 = huberWeight(res, gamma)
        # w0=w0/np.max(w0)
        w0 = w0.reshape(
            m,
        )
        D = np.diag(w0)
        adj = np.sum(w0)
        ts[k - 1, :] = thetaT  # storing all thetas
        ws[k - 1, :] = w0  # and weights
    if k == max_iter:  # if didn't converge
        ind = epsie.index(np.min(epsie))
        thetaT = ts[ind - 1, :]  # find most stable estimates
        w0 = ws[ind - 1, :]  # and weights
    return thetaT, w0


def gscreen(
    bx,
    by,
    SigmaUU,
    gamma,
    lamvec=np.linspace(0.01, 1.5, 15),
    rho=1,
    eps=1e-2,
    max_iter=15,
):
    if np.all(bx[:, 0] == 1):
        bx = bx[:, 1:]
        SigmaUU = SigmaUU[:, 1:][1:, :]
    # this function expects bx and by to already be transformed out of LD
    m = bx.shape[0]
    p = bx.shape[1]
    bics = []
    for ll in lamvec:
        thetaT, w0 = admmHuber(
            bx, by, SigmaUU, gamma=gamma, lam=ll, rho=rho, eps=eps, max_iter=max_iter
        )
        w0 = w0.squeeze()
        # w0=w0/np.max(w0)
        adj = np.sum(w0)
        res = by.squeeze() - bx @ (
            thetaT.reshape(
                p,
            )
        )
        rss = np.sum((w0 * res) ** 2)
        df2 = np.sum(thetaT != 0)
        if df2 == 0:
            break
        toadd = adj / m * np.log(rss / m) - np.log(adj) * m * df2
        # toadd=adj*np.log(rss)+(np.log(m-adj)+np.log(p))*df2
        bics.append(toadd)
    if len(bics) > 0:
        ind = np.argmin(np.array(bics).squeeze())
        thetaT, w0 = admmHuber(
            bx,
            by,
            SigmaUU,
            gamma=gamma,
            lam=lamvec[ind],
            rho=rho,
            eps=eps,
            max_iter=max_iter,
        )
        return thetaT, w0, bics[ind]
    else:
        return np.linspace(0, 0, p), np.linspace(1, 1, m), 0
