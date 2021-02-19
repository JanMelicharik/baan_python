import numpy as np
from math import gamma

def nu_lambda_post(nu_lam, nu_lam_0, lam):
    n = lam.shape[0]
    eta = 1 / nu_lam_0 + 0.5 * np.sum(lam - np.log(lam))
    return n * nu_lam / 2  * np.log(nu_lam) - n * np.log(gamma(nu_lam / 2)) - eta * nu_lam
