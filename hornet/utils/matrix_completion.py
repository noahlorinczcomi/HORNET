def ffimpute(X, K, dg, P, W, iter_max=25, eps=0.001):
    n, p = X.shape
    Del = numpy.isnan(X)
    # initialize
    X1 = X.copy()
    X1[Del] = 0
    S = P @ X1.T @ W @ X1 @ P
    Gamma = poet(S, K)
    Factor = X1 @ P @ Gamma
    L = Factor @ Gamma.T
    X1[Del] = L[Del]
    Gamma1 = Gamma * 0
    iter_ = 0
    while (numpy.linalg.norm(Gamma - Gamma1, ord="fro") > (eps * p**0.5)) & (
        iter_ < iter_max
    ):
        Gamma1 = Gamma.copy()
        S = P @ X1.T @ W @ X1 @ P
        Gamma = poet(S, K)
        Factor = X1 @ P @ Gamma
        L = Factor @ Gamma.T
        X1[Del] = L[Del]
        iter_ = iter_ + 1
    return L, Factor, Gamma


def ffimputeBIC(X, Kmax, dg, dim_proj, W, iter_max=15, eps=0.001, tuning=0.001):
    if dim_proj > 0:
        P = projectm(dg, k=dim_proj, tuning=tuning)
    if dim_proj == 0:
        P = makeI(X.shape[1])
    n, p = X.shape
    Del = numpy.isnan(X)
    Del = Del + 1 - 1  # converts to numeric (1=missing, 0=not)
    notDel = 1 - Del.copy()
    neff = numpy.sum(notDel)
    X1 = X.copy()
    X1[Del == 1] = 0
    bic = []
    for kk in range(1, Kmax + 1):
        L, Factor, Gamma = ffimpute(
            X=X, K=kk, P=P, W=W, dg=dg, iter_max=iter_max, eps=eps
        )
        df = n * kk + kk * p
        toAdd = (
            numpy.log(numpy.linalg.norm(X1 @ P - notDel * L, ord="fro") ** 2)
            + numpy.log(neff / (n + p)) * df / neff
        )
        bic.append(toAdd)
    kopt = bic.index(min(bic))
    L, Factor, Gamma = ffimpute(
        X=X, K=kopt, P=P, W=W, dg=dg, iter_max=iter_max, eps=eps
    )
    return L, Factor, Gamma
