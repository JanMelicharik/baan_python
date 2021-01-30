from numpy.random import gamma
from numpy import reshape

'''
Funkce pro generovani nahodnych cisel z G(mu, nu)
- podle Koop(2003)

    M ... pocet radku
    N ... pocet sloupcu
'''

def gamrnd_Koop(mu, nu, M, N):
    A = nu/2
    B = 2 * mu/nu

    # tato funkce vraci vektor o delce M x N
    y = gamma(A, B, M*N)

    Y = reshape(y,(M,N))

    return Y
