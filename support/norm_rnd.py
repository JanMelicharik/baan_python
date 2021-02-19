from scipy.linalg import cholesky as chol
from numpy.random import normal as randn
import numpy as np

'''
Nahodny vektor vice promennych dany matici sigma
    sigma ... ctvercova a symetricka kovariancni matice

Vystupem je nahodny vektor vyberu z normalniho rozlozeni
'''

def norm_rnd(sigma):
    h = chol(sigma)
    if isinstance(sigma, np.ndarray):
        rv = randn(size=(sigma.shape[0], 1))
    else:
        rv = randn(size=1)
    return  h.T @ rv
