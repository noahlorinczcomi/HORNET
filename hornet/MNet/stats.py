# standard imports
import numpy as np
import scipy as sp

# custom imports


def entropyloss(A, B, eps=1e-8):
    """
    Docstring for entropyloss
    
    :param A: Description
    :param B: Description
    :param eps: Description
    """
    C = A @ B
    a = np.trace(C)
    b, o = np.linalg.eig(C)
    b = np.real(b)
    mask = b > eps
    b = sum(np.log(b[mask]))
    return a - b - A.shape[1]


def Penamatrix(H2):
    """
    Docstring for Penamatrix
    
    :param H2: Description
    """
    S = np.eye(len(H2)) * 0
    for i in range(0, len(H2)):
        for j in range(0, len(H2)):
            S[i, j] = 1 / (H2[i] * H2[j]) ** 0.5
            S[j, i] = S[i, j]
    np.fill_diagonal(S, 0)
    return S


def spearmancov(A):
    p = A.shape[1]
    s = np.arange(0, p)
    # s=sp.stats.median_abs_deviation(A,axis=0)
    s = np.std(A, axis=0)
    s = np.array(s)
    R = sp.stats.spearmanr(A, axis=0).correlation
    R = 2 * np.sin(R * np.pi / 6)
    S = np.diag(s) @ R @ np.diag(s)
    return S


def soft(a, b):
    """
    Docstring for soft
    
    :param a: Description
    :param b: Description
    """
    c = abs(a) - b
    mask = c < 0
    c[mask] = 0
    return c * np.sign(a)


def mcp(a, lam, ga=3):
    """
    Docstring for mcp
    
    :param a: Description
    :param lam: Description
    :param ga: Description
    """
    b = abs(a)
    z = soft(a, lam) / (1 - 1 / ga)
    mask = b > (ga * lam)
    z[mask] = a[mask]
    return z


def MCPthreshold(S, lam, ga=3):
    """
    Docstring for MCPthreshold
    
    :param S: Description
    :param lam: Description
    :param ga: Description
    """
    p = S.shape[0]
    d = np.diag(S) ** 0.5
    R = cov2cor(S)
    s = R.reshape((p * p,))
    s1 = abs(s)
    s2 = mcp(s1, lam, ga)
    s2 = s2 * np.sign(s)
    S1 = s2.reshape((p, p))
    np.fill_diagonal(S1, 1)
    D = np.diag(d)
    S2 = D @ S1 @ D
    return S2


def cov2cor(A):
    """
    Docstring for cov2cor
    
    :param A: Description
    """
    A_ = A
    d = np.diag(A_)
    mask = d < 0
    d[mask] = 1
    D = np.diag(1 / d**0.5)
    return D @ A @ D


def matrixsqrt(A):
    """
    Docstring for matrixsqrt
    
    :param A: Description
    """
    d, v = np.linalg.eig(A)
    d1 = d * 0
    mask = d > 0
    d1[mask] = 1 / d[mask]
    d = d**0.5
    d1 = d1**0.5
    A = v @ np.diag(d) @ v.T
    B = v @ np.diag(d1) @ v.T
    return A, B


def positiveadj(A, min_eps=0.001):
    """
    Docstring for positiveadj
    
    :param A: Description
    :param min_eps: Description
    """
    d, v = np.linalg.eig(A)
    mask = d < min_eps
    d[mask] = 0
    B = v @ np.diag(d) @ v.T
    return B


def rowMedian(A):
    """
    Docstring for rowMedian
    
    :param A: Description
    """
    g = [np.median(A[_, :]) for _ in range(0, A.shape[0])]
    return g


def entropy_mcp_spearman_sampling(
    BETA,
    SE,
    Rnoise,
    lamvec=np.linspace(2, 15, 10) / 100,
    max_eps=0.001,
    max_iter=25,
    rho=0.25,
    min_eig=0.01,
    subtime=100,
    subfrac=0.66,
    subthres=0.95,
):
    """
    Docstring for entropy_mcp_spearman_sampling
    
    :param BETA: Description
    :param SE: Description
    :param Rnoise: Description
    :param lamvec: Description
    :param max_eps: Description
    :param max_iter: Description
    :param rho: Description
    :param min_eig: Description
    :param subtime: Description
    :param subfrac: Description
    :param subthres: Description
    """
    m, p = BETA.shape
    #########################tuning parameter selection##############################
    subvecerror = np.zeros((len(lamvec), subtime))
    Thetalist = np.zeros((len(lamvec), subtime, p, p))

    for j in range(0, subtime):
        indsub = np.random.choice(m, size=round(subfrac * m), replace=False)
        mask = np.zeros((m,), dtype="bool")
        mask[indsub] = True
        BETAS = BETA[mask]
        SES = SE[indsub]
        S = spearmancov(BETAS) * BETAS.shape[0]
        M = S * 0
        for ii in range(0, BETAS.shape[0]):
            M = M + SES[ii, :] * (Rnoise * SES[ii, :]).T
        S1 = cov2cor(S - M)
        Theta0 = np.linalg.inv(MCPthreshold(S1, 0.1))
        BETAS2 = BETA[mask == False, :]
        SES2 = SE[mask == False, :]
        S2 = spearmancov(BETAS2) * BETAS2.shape[0]
        M2 = S2 * 0
        for ii in range(0, BETAS2.shape[0]):
            M2 = M2 + SES[ii, :] * (Rnoise * SES2[ii, :]).T
        S2 = cov2cor(S2 - M2)

        for i in range(0, len(lamvec)):
            Theta = Theta0.copy()
            Theta1 = Theta0 * 0
            Delta1 = Theta.copy()
            Gamma1 = Delta1 * 0
            Delta2 = Theta.copy()
            Gamma2 = Delta2 * 0
            error = np.linalg.norm(Theta - Theta1, ord="fro") / p**0.5
            iter = 0
            while (error > max_eps) & (iter < max_iter):
                Theta1 = Theta.copy()
                Q = S1 + Gamma1 - rho * Delta1
                Theta = (-Q + matrixsqrt(Q @ Q + 4 * rho * np.eye(p))[0]) / 2 / rho
                Delta1 = mcp(
                    (Theta + Gamma1 / rho).reshape((p * p,)), lamvec[i] / rho, ga=3
                ).reshape((p, p))
                Gamma1 = Gamma1 + rho * (Theta - Delta1)
                if (iter % 5 == 0) & (min(np.linalg.eig(Theta)[0]) < min_eig):
                    Q = S1 + Gamma1 + Gamma2 - rho * (Delta1 + Delta2)
                    Theta = (
                        (-Q + matrixsqrt(Q @ Q + 8 * rho * np.eye(p))[0]) / 4 / rho
                    )
                    Delta1 = mcp(
                        (Theta + Gamma1 / rho).reshape((p * p,)), lamvec[i] / rho, ga=3
                    ).reshape((p, p))
                    Gamma1 = Gamma1 + rho * (Theta - Delta1)
                    Delta2 = positiveadj(Theta + Gamma2 / rho, min_eps=min_eig)
                    Gamma2 = Gamma2 + rho * (Theta - Delta2)
                # end if
                iter = iter + 1
                error = np.linalg.norm(Theta - Theta1, ord="fro") / p**0.5
            # end while
            df = sum((sum(Delta1 != 0) + p) / 2)
            subvecerror[i, j] = entropyloss(S2, Delta1) + np.log(m) / m * df
            Theta0 = Delta1.copy()
            Thetalist[i, j, :, :] = Theta0
        # end for lamvec
    # end for subtime
    istar = np.argmin(rowMedian(subvecerror))
    Thetalist = Thetalist[istar, :, :, :]
    K = Thetalist[1, :, :] * 0
    Thetasub = K.copy()
    for i in range(0, subtime):
        K = K + sum(sum(Thetalist[i, :, :] != 0)) / subtime
        Thetasub = Thetasub + Thetalist[i, :, :] / subtime
    S = spearmancov(BETA) * BETA.shape[0]
    M = S * 0
    for i in range(0, BETA.shape[0]):
        M = M + SE[i, :] * (Rnoise * SE[i, :]).T
    S1 = cov2cor(S - M)
    Theta0 = np.linalg.inv(MCPthreshold(S1, 0.1))
    Theta = Theta0.copy()
    Theta1 = Theta0 * 0
    Delta1 = Theta.copy()
    Gamma1 = Delta1 * 0
    Delta2 = Theta.copy()
    Gamma2 = Delta2 * 0
    error = np.linalg.norm(Theta - Theta1, ord="fro") / p**0.5
    iter = 0
    while (error > max_eps) & (iter < (2 * max_iter)):
        Theta1 = Theta.copy()
        Q = S1 + Gamma1 - rho * Delta1
        Theta = (-Q + matrixsqrt(Q @ Q + 4 * rho * np.eye(p))[0]) / 2 / rho
        Delta1 = mcp(
            (Theta + Gamma1 / rho).reshape((p * p,)), lamvec[istar] / rho, ga=3
        ).reshape((p, p))
        Gamma1 = Gamma1 + rho * (Theta - Delta1)
        if (iter % 5 == 0) & (min(np.linalg.eig(Theta)[0]) < min_eig):
            Q = S1 + Gamma1 + Gamma2 - rho * (Delta1 + Delta2)
            Theta = (-Q + matrixsqrt(Q @ Q + 8 * rho * np.eye(p))[0]) / 4 / rho
            Delta1 = mcp(
                (Theta + Gamma1 / rho).reshape((p * p,)), lamvec[istar] / rho, ga=3
            ).reshape((p, p))
            Gamma1 = Gamma1 + rho * (Theta - Delta1) * (K > subthres)
            Delta2 = positiveadj(Theta + Gamma2 / rho, min_eps=min_eig)
            Gamma2 = Gamma2 + rho * (Theta - Delta2)
        # end for
        iter = iter + 1
        error = np.linalg.norm(Theta - Theta1, ord="fro") / p**0.5
    # end while
    Theta0 = Delta1.copy()
    A = dict()
    A["Theta"] = Theta0
    A["Theta_Alt"] = Thetasub
    A["K"] = K
    A["cv_error"] = subvecerror
    return A
