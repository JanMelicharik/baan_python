import numpy as np

from math import sqrt, log, pi
from numpy.linalg import inv, det
from numpy import transpose as t, diag
from scipy.special import gammaln

'''
Funkce pro bayesianskou analyzu NLRM (prirozene konjugovana apriorni hustota)

Podminky: y, x, beta_0 a v_0 musi byt matice

Vstupy:
    y ........................ vektor vysvetlujici promenne velikosti N x 1
    x ........................ matice planu velikosti N x k
    beta_0, v_0, h_0, nu_0 ... apriorni parametry
                               p(beta, h) ~ NG(beta_0, v_0, h_0, nu_0)
    beta_1, v_1, h_1, nu_1 ... aposteriorni hyperparametry
                               p(beta, h|y) ~ NG(beta_0, v_0, h_0, nu_0)                 

Vystupy:
    b0_cov ........... apriorni kovariancni matice
    b0_std, h0_std ... vektory apriornich smerodatnych odchylek parametru
    b1_cov ........... posteriorni kovariancni matice
    b1_std, h1_std ... vektory posteriornich smerodatnych odchylek parametru
    log_ML ........... logaritmus marginalni verohodnosti modelu
'''

def my_nlrm(y, x, beta_0, v_0, h_0, nu_0):
    b0_cov = nu_0 * h_0**(-1) / (nu_0 - 2) * v_0                # analogie (3.16) z Koop (2003)
    b0_std = np.array([[sqrt(var)] for var in diag(b0_cov)])
    h0_std = sqrt(2 * h_0 / nu_0)                               # analogie odmocniny z (3.19) z Koop (2003)

    n = len(y)      # pocet pozorovani

    b_ols  = inv(t(x) @ x) @ t(x) @ y                           # (3.5) z Koop (2003)
    nu_ols = n - x.shape[1]                                     # (3.4) z Koop (2003)
    s2_ols = 1 / nu_ols * t(y - x @ b_ols) @ (y - x @ b_ols)    # (3.6) z Koop (2003)

    v_1    = inv(inv(v_0) + t(x) @ x)                           # (3.9) z Koop (2003)
    beta_1 = v_1 @ (inv(v_0) @ beta_0 + t(x) @ x @ b_ols)       # (3.11) z Koop (2003)
    nu_1   = nu_0 + n                                           # (3.12) z Koop (2003)

    h_1 = nu_1 / (nu_0 / h_0 + nu_ols * s2_ols + t(b_ols - beta_0) @ inv(v_0 + inv(t(x) @ x)) @ (b_ols - beta_0))
    # (3.13) z Koop (2003), kdy h_1 = s2_1^-1

    b1_cov = nu_1 / h_1 / (nu_1 - 2) * v_1                      # (3.16) z Koop (2003)
    b1_std = np.array([[sqrt(var)] for var in diag(b1_cov)])
    h1_std = sqrt(2 * h_1 / nu_1)                               # odmocnina z (3.19) z Koop (2003)

    log_c  = gammaln(nu_1 / 2) + nu_0 / 2 * log(nu_0 / h_0) - gammaln(nu_0 / 2) - n / 2 * log(pi)
    # logaritmus (3.35) z Koop (2003)
    log_ml = log_c + 1/2 * (log(det(v_1)) - log(det(v_0))) - nu_1 / 2 * log(nu_1 / h_1)
    # logaritmus (3.34) z Koop (2003)

    results = {
        'x'      : x,
        'y'      : y,
        'beta_0' : beta_0,
        'v_0'    : v_0,
        'h_0'    : h_0,
        'nu_0'   : nu_0,
        'b0_cov' : b0_cov,
        'b0_std' : b0_std,
        'h0_std' : h0_std,
        'n'      : n,
        'b_ols'  : b_ols,
        'nu_ols' : nu_ols,
        's2_ols' : s2_ols,
        'beta_1' : beta_1,
        'h_1'    : h_1,
        'v_1'    : v_1,
        'nu_1'   : nu_1,
        'b1_cov' : b1_cov,
        'b1_std' : b1_std,
        'h1_std' : h1_std,
        'log_c'  : log_c,
        'log_ml' : log_ml
    }

    return results
