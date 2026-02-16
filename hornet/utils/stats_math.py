def projectm(dg, bdeg=3, k=10, tuning=0):
    q1 = numpy.quantile(dg, 0.05)
    x = dg - numpy.min(dg)
    x = x / numpy.max(x)
    B, D = cosinebase(x, k)
    P = B @ numpy.linalg.inv(B.T @ B + tuning * D / p) @ B.T
    return P


def cosinebase(t, k):
    p = len(t)
    A = numpy.ones((p, k + 1))
    d = [0]
    for _ in range(1, k + 1):
        # A=numpy.column_stack((A,numpy.sqrt(2)*numpy.cos(numpy.pi*_*t)))
        A[:, _] = numpy.sqrt(2) * numpy.cos(numpy.pi * _ * t)
        d.append(_**4)
    D = numpy.zeros((k + 1, k + 1))
    numpy.fill_diagonal(D, d)
    return A, D


def CD(
    x,
    y,
    UU,
    UV,
    VV,
    Rinv,
    lambda_,
    alpha_,
    eps=0.00001,
    max_iter=100,
    norm_stoppage=100,
    pen="scad",
):
    p = x.shape[1]
    m = x.shape[0]
    I = makeI(p)
    # no adjustment for horizontal pleiotropy here
    est0 = numpy.linalg.inv(x.T @ Rinv @ x + lambda_ * I - m * UU) @ (
        x.T @ Rinv @ y - m * UV
    )
    nns = []
    nn = 100
    for cycle_ in range(0, max_iter):
        normPrev = nn
        bjhat = []
        for p_ in range(0, p):
            xp_ = x[:, p_]
            r_j = y - numpy.delete(x, obj=p_, axis=1) @ numpy.delete(
                est0, obj=p_, axis=0
            )
            rx = (1 / m) * (r_j.T @ xp_)
            rx = rx[0]  # turn to int
            dx = (1 / m) * (xp_.T @ xp_) + (1 - alpha_) * lambda_
            if pen.lower() == "lasso":
                toAdd = Sop(rx, lambda_) / dx
            else:
                toAdd = scad(numpy.array(rx), 3.7, lambda_) / dx
            bjhat.append(toAdd)
        bjhat = numpy.array(bjhat)
        bjhat = bjhat.reshape((bjhat.shape[0], 1))
        nn = ((est0 - bjhat).T @ (est0 - bjhat)) ** 2
        nn = nn[0][0]
        nns.append(nn)
        est0 = bjhat
        # if norm(est0-bjhat)>norm_stoppage and this norm>previous norm, just stop - not going to converge
        if nn < (eps * p) or ((nn > norm_stoppage and nn > normPrev)):
            break
    # end loops
    dd = {"est": est0, "k": cycle_}
    return dd
