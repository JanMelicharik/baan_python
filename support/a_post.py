import numpy as np
from numpy.linalg import inv
from math import log

def a_post(a, b, h, y, x, z):
    # n - pocet pozorovani; y musi byt sloupcova matice nebo vektor
    n = y.shape[0]
    diagonal = ((z @ a + 1) ** 2).ravel()
    omega = np.diag(diagonal)
    pom = np.sum([log(val) for val in diagonal])
    f = y - x @ b
    return - 0.5 * pom - 0.5 * h * f.T @ inv(omega) @ f
