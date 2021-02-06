# Nastaveni cesty do domovske slozky
import sys
import pdb
sys.path.insert(1, "/".join(sys.path[0].split("/")[:-1]))

# importy pomocnych funkci
from support.progress_info import progress_bar
from support.gamm_rnd_koop2 import gamm_rnd_koop2
from support.norm_rnd import norm_rnd
from support.log_post_ces import log_post_ces
from support.lik_ces import lik_ces
from support.prior_ces import prior_ces

from numpy.random import multivariate_normal as mvn
from numpy.linalg import (inv, det)
from numpy.random import uniform
from numpy.random import normal
from scipy.stats.distributions import chi2
from math import (log, pi, exp)
from pandas import read_csv
from tabulate import tabulate

import warnings
import math
import numpy as np
import matplotlib.pyplot as plt

# V prubehu vypoctu mohou nastat dva typy warningu - nemaji vliv na vysledek scriptu
warnings.filterwarnings("ignore", r"overflow encountered in (power|matmul)")

data = read_csv("mexico.csv", delimiter=",")

# ===== 1. Priprava dat (do podoby indexu) =====

y = np.array([
        data["gdp"] / np.mean(data["gdp"])
    ]).T

x = np.array([
        [1]* len(y),
        data["labor"] / np.mean(data["labor"]),
        data["capital"] / np.mean(data["capital"]),
    ]).T

n = y.shape[0]

# ===== 2. Apriorni hustoty a apriorni hyperparametry =====
# p(gamma) ~ N(gamma_0, v_0)
gamma_0 = np.array(
    [
        [1],
        [0.5],
        [0.5],
        [1]
    ]
)
k = gamma_0.size       # Pocet parametru modelu
v_0 = np.diag(
    [
        0.5**2,
        0.25**2,
        0.25**2,
        0.5**2
    ]
)
# p(h) ~ G(h_0, nu_0)
h_0 = 1 / 0.5**2
nu_0 = 5

# ===== 3. Metropolis within Gibbs - nastaveni =====
s = 50_000 + 1
s_0 = 30_000 + 1
s_1 = s - s_0

# Ukladani vzorku
gamma = np.zeros((k, s))
gamma[:,[0]] = gamma_0
h = np.array([0.0] * s)

# Nastaveni RW M-H algoritmu:
# Kandidatska hustota ~ N(gamma(s-1), sigma)
d = 0.5                         # Skalovaci konstanta
# sigma = d * np.identity(k)      # 
sigma = d * np.array(
    [
        [ 0.0680, -0.0343, -0.0284, -0.0024],
        [-0.0343,  0.0449, -0.0021,  0.0037],
        [-0.0284, -0.0021,  0.0341, -0.0015],
        [-0.0024,  0.0037, -0.0015,  0.2144],
    ]
)
count = 0                       # Pocitadlo akceptovasnych vzorku
nu_1 = nu_0 + n                 # 5.23
log_A = h = np.array([0.0] * s)
# ===== 4. Metropolis within Gibbs =====
print("Metropolis within Gibbs ...")
for i in range(1,s):
    # a) podminena hustota p(h|gamma,y) ~ G(h_1,nu_1)
    f_x = gamma[0][i-1] * \
        ( gamma[1][i-1] * x[:,[1]]**gamma[3][i-1] + gamma[2][i-1] * x[:,[2]]**gamma[3][i-1] ) ** \
        ( 1 / gamma[3][i-1] )
    h_1 = (1 / nu_1 * ((y - f_x).T @ (y - f_x) + nu_0 / h_0))**(-1)
    h[i] = gamm_rnd_koop2(h_1, nu_1)
    # b) podminena hustota p(gamma|h,y) -> logaritmus jadrove podminene hustoty = externi funkce log_post_ces.py
    # Generovani kandidatu z kandidatske hustoty
    gamma_ast = gamma[:,[i-1]] + norm_rnd(sigma)
    log_accept = min(
        log_post_ces(y, x, gamma_ast, h[i], gamma_0, v_0) - \
        log_post_ces(y, x, gamma[:,[i-1]], h[i], gamma_0, v_0),
        0
    )
    # Rozhodnuti o akceptaci
    if log_accept > log(uniform()):
        gamma[:,[i]] = gamma_ast
        count += 1
    else:
        gamma[:,[i]] = gamma[:,[i-1]]

    progress_bar(i,s)

# ===== 5. Prezentace vysledku a posteriorni analyza =====

# a) Vyhozeni prvnich s_0 vzorku
gamma = gamma[:,s_0:]
h = h[s_0:]

# b) Vypocet aposteriornich momentu
e_gamma = np.mean(gamma, axis=1)
d_gamma = np.var(gamma, axis=1)
std_gamma = np.std(gamma, axis=1)
e_h = np.mean(h)
d_h = np.var(h)
std_h = np.std(h)

# Apriorni momenty
e_gamma_0 = gamma_0
d_gamma_0 = np.diag(v_0)
std_gamma_0 = [math.sqrt(var) for var in d_gamma_0]
e_h_0 = h_0
d_h_0 = 2 * e_h_0**2 / nu_0
std_h_0 = math.sqrt(d_h_0)

headers = ["Parametr", "E_prior", "std_prior", "E_post", "std_post"]
table = []
for i in range(len(e_gamma)):
    table.append([
        f"gamma_{i+1}",
        round(e_gamma_0[i][0], 4),
        round(std_gamma_0[i], 4),
        round(e_gamma[i], 4),
        round(std_gamma[i], 4)])
table.append(["h", round(e_h_0, 4), round(std_h_0, 4), round(e_h, 4), round(std_h, 4)])

print("\nApriorni a aposteriorni parametry:")
print(tabulate(table, headers, tablefmt="pretty"))

# ===== 6. Vypocet marginalni verohodnosti modelu metodou Gelfanda a Deye =====
# a) Zjednodusena varianta
gd_simple = 0
print("Statistika Gelfanda a Deye (zjednodusene) ...")
for i in range(s_1):
    gd_simple = gd_simple + 1 / s_1 * 1 / lik_ces(y, x, gamma[:,[i]], h[i])
    progress_bar(i,s_1)

print("\nMarginalni verohodnost CES produkcni funkce (zjednodusena varianta):")
print(f"Marginal likelihood = {1/gd_simple}")
print(f"Marginal likelihood (log) = {math.log(1/gd_simple)}\n")

# b) plna varianta
# vektor 'h' je jednorozmerny - aby bylo mozne jej pridat jako radkovou matici pod
# matici gamma, je adresovan pomoci druhe dimenze (None) - viz. NumPy, Broadcasting
theta = np.r_[gamma, h[None,:]]
theta_hat = np.mean(theta, axis=1)
sigma_hat = np.cov(theta)
kk = theta_hat.size     # protoze jde o vektor, lze pouzit size namisto shape
pp = 0.01               # Pro (1 - p) procentni kvantil chi-kvadrat rozdeleni
chi_pp = chi2.ppf(1 - pp, kk)

# Simulace integracni konstanty pro omezenou apriorni hustotu
count_g = 0
for i in range(s):
    pom = gamma_0 + norm_rnd(v_0)   # p(gamma) ~ N(gamma_0,V_0)
    count_g += min(pom) > 0         

fth = np.array([0.0] * s)
pri = np.array([0.0] * s)
lik = np.array([0.0] * s)

int_c = s / count_g
gd = 0      # Statistika Gelfanda a Deye
print("Statistika Gelfanda a Deye ...")
for i in range(s_1):
    t_theta = (theta[:,[i]] - theta_hat[:, None]).T @ inv(sigma_hat) @ (theta[:,[i]] - theta_hat[:, None])
    # Funkce hustoty vicerozmerneho rozdeleni
    f_theta = 1 / (1 - pp) * 1 / (2 * pi) ** (kk / 2) * \
            det(sigma_hat) ** (- 1 / 2) * exp(- 1 / 2 * t_theta) \
            if t_theta <= chi_pp else 0
    prior = prior_ces(gamma[:,[i]], h[i], gamma_0, v_0, h_0, nu_0)
    like = lik_ces(y, x, gamma[:,[i]], h[i])
    lik[i] = like
    pri[i] = prior
    fth[i] = f_theta
    gd += f_theta / (int_c * prior * like) * 1 / s_1

    progress_bar(i,s_1)

print("\nMarginalni verohodnost CES produkcni funkce (plna varianta):")
print(f"Marginal likelihood = {1/gd}")
print(f"Marginal likelihood (log) = {math.log(1/gd)}\n")

# # ===== 7. Srovnani skutecnych a modelovych momentu =====
e_y_ast = np.zeros((1, s_1))    # Vektor modelovych strednich hodnot
d_y_ast = np.zeros((1, s_1))    # Vektor modelovych rozptylu
std_y_ast = np.zeros((1, s_1))  # Vektor modelovych sm. odchylek

# # Simulace (generovani umelych dat)
for i in range(s_1):
    f_xx = gamma[0][i] * \
        ( gamma[1][i] * x[:,[1]] ** gamma[3][i] + gamma[2][i] * x[:,[2]] ** gamma[3][i] ) ** \
        ( 1 / gamma[3][i] )
    y_ast = f_xx + normal() * math.sqrt(1 / h[i])
    e_y_ast[0, [i]] = np.mean(y_ast)
    d_y_ast[0, [i]] = np.var(y_ast)
    std_y_ast[0, [i]] = np.std(y_ast)

a = np.sum(e_y_ast < np.mean(y))
b = np.sum(d_y_ast < np.var(y))
c = np.sum(std_y_ast < np.std(y))

# Vypocet predikcnich (jednostrannych) p-hodnot
# a) pro stredni hodnotu
p_e = (a / s_1) * ((a / s_1) <= 0.5) + (1 - a / s_1) * ((a / s_1) > 0.5)
# b) pro rozptyl
p_d = (b / s_1) * ((b / s_1) <= 0.5) + (1 - b / s_1) * ((b / s_1) > 0.5)
# c) pro sm. odchylku
p_std = (c / s_1) * ((c / s_1) <= 0.5) + (1 - c / s_1) * ((c / s_1) > 0.5)

print("Predikcni p-hodnoty:")
print(f" p_E = {round(p_e, 4)}\n p_D = {round(p_d, 4)}\n p_std = {round(p_std, 4)}")

fix, ax = plt.subplots(ncols=3, figsize=(20,5))
ax[0].set_title("Simulovane stredni hodnoty")
ax[0].hist(e_y_ast[0], bins=30, histtype="bar", rwidth=0.9)
ax[0].axvline(np.mean(y), color="r", linewidth=2, linestyle="dashed")
ax[1].set_title("Simulovane rozptyly")
ax[1].hist(d_y_ast[0], bins=30, histtype="bar", rwidth=0.9)
ax[1].axvline(np.var(y), color="r", linewidth=2, linestyle="dashed")
ax[2].set_title("Simulovane sm. odchylky")
ax[2].hist(std_y_ast[0], bins=30, histtype="bar", rwidth=0.9)
ax[2].axvline(np.std(y), color="r", linewidth=2, linestyle="dashed")

plt.show(block=False)
pdb.set_trace()
