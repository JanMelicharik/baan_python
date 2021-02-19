from scipy.linalg import cholesky as chol
from numpy.random import normal as randn
import numpy as np
from numpy.linalg import eig
import sys
from termcolor import colored

def wish_rnd(sigma: np.ndarray, v: int):
    n, k = sigma.shape
    eigenvals = eig(sigma)
    if not all([val > 0 for val in eigenvals]):
        print(colored("wish_rnd: Matice musi byt positivne definitivni.", "red"))
        sys.exit(1)

    t = chol(sigma)
    y = t.T @ randn(n, v)
    return y @ y.T

    