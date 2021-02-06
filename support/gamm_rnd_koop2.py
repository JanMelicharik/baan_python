from numpy.random import gamma

def gamm_rnd_koop2(mu, nu, m=1):
    a = nu / 2
    b = 2 * mu / nu
    return b * gamma(a, size=m)
    
