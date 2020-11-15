from scipy.linalg import cholesky as chol
from numpy import shape
from numpy import transpose as t
from numpy.random import normal as randn

'''
Nahodny vektor vice promennych dany matici sigma
    sigma ... ctvercova a symetricka kovariancni matice

Vystupem je nahodny vektor vyberu z normalniho rozlozeni
'''

def norm_rnd(sigma):
    h = chol(sigma)
    size = shape(sigma)
    rv = randn(size = (size[0],1))
    y = t(h) @ rv

    return y



