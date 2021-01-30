from pandas import read_excel, DataFrame
from math import sqrt, lgamma, log, pi
from numpy import diag, ones, dot, array, asmatrix, shape
from numpy import transpose as t
from numpy.linalg import inv, det

'''
Funkce pro bayesianskou analyzu NLRM (prirozene konjugovana apriorni hustota)

Vstupy:
    y ........................ vektor vysvetlujici promenne velikosti N x 1
    X ........................ matice planu velikosti N x k
    beta_0, V_0, h_0, nu_0 ... apriorni parametry
                               p(beta, h) ~ NG(beta_0, V_0, h_0, nu_0)
    beta_1, V_1, h_1, nu_1 ... aposteriorni hyperparametry
                               p(beta, h|y) ~ NG(beta_0, V_0, h_0, nu_0)                 

Vystupy:
    b0_cov ........... apriorni kovariancni matice
    b0_std, h0_std ... vektory apriornich smerodatnych odchylek parametru
    b1_cov ........... posteriorni kovariancni matice
    b1_std, h1_std ... vektory posteriornich smerodatnych odchylek parametru
    log_ML ........... logaritmus marginalni verohodnosti modelu
'''

def my_nlrm(y, X, beta_0, V_0, h_0, nu_0):
# Vypocet charakteristik apriornich hustot
    # Vypocet kovariancni matice pro beta
    b0_cov = nu_0 * (1/h_0)/(nu_0 - 2) * V_0    # analogie (3.16) z Koop (2003)
    b0_std = [sqrt(cov) for cov in diag(b0_cov)]

    # Apriorni sm. odchylka pro h
    h0_std = sqrt(2 * h_0/nu_0)                 # analogie (3.19) z Koop (2003)
    
    # Vypocet aposteriornich hyperparametru
    N = len(y)

    # Odhady OLS
    b_OLS = inv(t(X) @ X) @ t(X) @ y                            # (3.5) z Koop (2003)
    nu_OLS = N - shape(X)[1]                                    # (3.4) z Koop (2003)
    s2_OLS = t(y - dot(X, b_OLS)) @ (y - dot(X, b_OLS))/nu_OLS  # (3.6) z Koop (2003)

    # Aposteriorni hyperparametry
    V_1 = inv(inv(V_0) + (t(X) @ X))                            # (3.10) z Koop (2003)
    beta_1 = V_1 @ inv(V_0) @ beta_0 + t(X) @ dot(X, b_OLS)     # (3.11) z Koop (2003)
    beta_1 = asmatrix(beta_1)
    nu_1 = nu_0 + N                                             # (3.12) z Koop (2003)
    h_1 = nu_1 * 1/(nu_0 * 1/h_0 + nu_OLS * s2_OLS + \
          t(b_OLS - beta_0) @ inv(V_0 + \
          inv(t(X) @ X)) @ (b_OLS - beta_0))                    # (3.13) z Koop (2003)

    # Aposteriorni kovariancni matice a smerodatne odchylky
    b1_cov = (nu_1 * 1/h_1)/(nu_1 - 2) * V_1                    # (3.16) z Koop (2003)
    b1_std = [sqrt(cov) for cov in diag(b1_cov)]

    # Aposteriorni smerodatna odchylka pro presnost chyby h
    h1_std = sqrt(2 * h_1/nu_1)                                 # (3.19) z Koop (2003), odmocneno

    # Logaritmus marginalni verohodnosti
    log_c = lgamma(nu_1/2) + nu_0/2 * log(nu_0/h_0) - \
        lgamma(nu_0/2) - N/2 * log(pi)                          # (3.35) z Koop (2003), logaritmovano
    log_ML = log_c + 1/2 * (log(det(V_1)) - log(det(V_0))) - \
             nu_1/2 * log(nu_1/h_1)                             # (3.34) z Koop (2003), logaritmovano

    # Ulozeni vysledku
    results = {
        'b0_cov' : b0_cov,
        'b0_std' : b0_std,
        'h0_std' : h0_std,
        'beta_1' : beta_1,
        'h_1'    : h_1,
        'V_1'    : V_1,
        'nu_1'   : nu_1,
        'b1_cov' : b1_cov,
        'b1_std' : b1_std,
        'h1_std' : h1_std,
        'log_ML' : log_ML
    }

    return results
