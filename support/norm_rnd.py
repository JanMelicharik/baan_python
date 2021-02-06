from scipy.linalg import cholesky as chol
from numpy.random import normal as randn

'''
Nahodny vektor vice promennych dany matici sigma
    sigma ... ctvercova a symetricka kovariancni matice

Vystupem je nahodny vektor vyberu z normalniho rozlozeni
'''

def norm_rnd(sigma):
    h = chol(sigma)
    rv = randn(size=(sigma.shape[0], 1))
    return  h.T @ rv
