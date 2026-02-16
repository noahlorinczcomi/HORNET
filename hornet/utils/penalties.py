# standard imports
from hornet.utils.stats import fPos
import numpy as np

# custom imports
from MNet import *  # this is a source file in the HORNET directory


def scad(z, a, alambda):
    cond1 = abs(z) <= (2 * alambda)
    cond2 = ((2 * alambda) < abs(z)) & (abs(z) <= (a * alambda))
    cond3 = abs(z) > (a * alambda)
    z[cond1] = np.sign(z[cond1]) * fPos(abs(z[cond1]) - alambda)
    z[cond2] = ((a - 1) * z[cond2] - np.sign(z[cond2]) * a * alambda) / (a - 2)
    z[cond3] = z[cond3]
    return z


def dSCAD(a, lam, gamma=3.7):
    a = abs(a)
    a = a.reshape((a.shape[0],))
    z = a.copy()
    mask1 = a < lam
    mask2 = a > lam
    mask3 = a > (gamma * lam)
    z[mask1] = lam
    z[mask2] = (gamma * lam - z[mask2]) / (gamma - 1)
    z[mask3] = 0
    return z


def soft(a, b):
    c = abs(a) - b
    mask = c < 0
    c[mask] = 0
    return c * np.sign(a)


def Sop(theta, lambda_):
    toReturn = []
    for jj in range(0, len(theta)):
        toS = np.abs(theta[jj]) - lambda_
        z = toS if toS > 0 else 0
        toReturn.append(np.sign(theta[jj]) * z)
    return toReturn


def mcp(a, lam, ga=3):
    b = abs(a)
    z = soft(a, lam) / (1 - 1 / ga)
    mask = b > (ga * lam)
    z[mask] = a[mask]
    return z


def softThresholdLD(R, lambda_):
    Rsoft = soft(R, lambda_)  # will NOT retain 1s on diagonal
    np.fill_diagonal(Rsoft, 1)  # retain correlation matrix structure
    return Rsoft


def huberWeight(e, gamma):
    e_ = np.array(e)
    mask = abs(e_) <= gamma
    w = e_ * 0
    w[mask] = 1
    w[mask == False] = gamma / abs(e[mask == False])
    return w
