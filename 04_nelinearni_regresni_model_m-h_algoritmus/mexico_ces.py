# Nastaveni cesty do domovske slozky
import sys
root_dir = "/".join(sys.argv[0].split("/")[:-3])
sys.path.insert(1, root_dir + "/Support/Python")

# importy pomocnych funkci
from support.progress_info import progress_bar
from support.gamm_rnd_koop import gamm_rnd_koop
from support.norm_rnd import norm_rnd
from support.log_post_ces import log_post_ces
from support.lik_ces import lik_ces
from support.prior_ces import prior_ces

from numpy import transpose as t
from numpy.linalg import (inv, det)
from numpy.random import uniform
from scipy.stats.distributions import chi2
from math import (log, pi, exp)
from pandas import read_csv
from tabulate import tabulate

import math
import numpy as np
import matplotlib.pyplot as plt

data = read_csv("mexico.csv", delimiter=",")

# ===== 1. Priprava dat (do podoby indexu) =====

y = np.array([
        data["gdp"] / np.mean(data["gdp"])
    ]).transpose()

x = np.array([
        [1]* len(y),
        data["labor"] / np.mean(data["labor"]),
        data["capital"] / np.mean(data["capital"]),
    ]).transpose()

n = y.shape[0]

# ===== 2. Apriorni hustoty a apriorni hyperparametry =====
# p(gamma) ~ N(gamma_0, v_0)
gamma_0 = np.array([[1],[0.5],[0.5],[1]])
k = gamma_0.shape[0]       # Pocet parametru modelu
v_0 = np.diag([0.5**2, 0.25**2, 0.25**2, 0.5**2])
# p(h) ~ G(h_0, nu_0)
h_0 = 1 / 0.5**2
nu_0 = 5

# ===== 3. Metropolis within Gibbs - nastaveni =====
s = 5_000 + 1
s_0 = 3_000 + 1
s_1 = s - s_0

# Ukladani vzorku
gamma = np.zeros((k, s))
gamma[:,0] = gamma_0[:,0]
h = [0] * s

# Nastaveni RW M-H algoritmu:
# Kandidatska hustota ~ N(gamma(s-1), sigma)
d = 0.5                         # Skalovaci konstanta
sigma = d * np.identity(k)      # Prvotni nastaveni kov. matice M-H algoritmu
count = 0                       # Pocitadlo akceptovanych vzorku
# sigma = d * np.array([
#     [ 0.0680, -0.0343, -0.0284, -0.0024],
#     [-0.0343,  0.0449, -0.0021,  0.0037],
#     [-0.0284, -0.0021,  0.0341, -0.0015],
#     [-0.0024,  0.0037, -0.0015,  0.2144]])

# ===== 4. Metropolis within Gibbs =====

nu_1 = nu_0 + n     # 5.23
for i in range(1,s):
# a) podminena hustota p(h|gamma,y) ~ G(h_1,nu_1)
    f_x = gamma[0][i-1] * (gamma[1][i-1] * x[:,1]**gamma[3][i-1] + gamma[2][i-1] * x[:,2]**gamma[3][i-1])**(1 / gamma[3][i-1])
    f_x = np.array([f_x]).transpose()   # Predchozi vypocet vrati vektor ve spatnem tvaru
    h_1 = 1 / (1 / nu_1 * (t(y - f_x) @ (y - f_x) + nu_0 * 1 / h_0))    # 5.23
    h_1 = h_1[0][0]
    h[i] = gamm_rnd_koop(h_1, nu_1, 1)[0]

    # b) podminena hustota p(gamma|h,y) -> logaritmus jadrove podminene hustoty = externi funkce log_post_ces.py
    # Generovani kandidatu z kandidatske hustoty
    gamma_ast = gamma[:,[i-1]] + norm_rnd(sigma)
    gamma_back = gamma[:,[i-1]]
    log_accept = min(log_post_ces(y, x, gamma_ast, h[i], gamma_0, v_0),
                    log_post_ces(y, x, gamma_back, h[i], gamma_0, v_0),
                    0)

    random_number = uniform()

    if log_accept > log(random_number):
        gamma[:,[i]] = gamma_ast
        count += 1
    else:
        gamma[:,[i]] = gamma[:,[i - 1]]

    progress_bar(i,s)

# ===== 5. Prezentace vysledku a posteriorni analyza =====

# a) Vyhozeni prvnich s_0 vzorku
gamma = np.delete(gamma, range(s_0), axis=1)
h = h[s_0:]

# b) Vypocet aposteriornich momentu
e_gamma = np.mean(gamma, axis=1)
d_gamma = np.var(gamma, axis=1)
sd_gamma = [math.sqrt(var) for var in d_gamma]
e_h = np.mean(h)
d_h = np.var(h)
sd_h = math.sqrt(d_h)

# Apriorni momenty
e_gamma_0 = gamma_0
d_gamma_0 = np.diag(v_0)
sd_gamma_0 = [math.sqrt(var) for var in d_gamma_0]
e_h_0 = h_0
d_h_0 = 2 * e_h_0**2 / nu_0
sd_h_0 = math.sqrt(d_h_0)

headers = ["Parametr", "E_prior", "std_prior", "E_post", "std_post"]
table = []
for i in range(len(e_gamma)):
    table.append([
        f"gamma_{i+1}",
        round(e_gamma[i], 4),
        round(sd_gamma[i], 4),
        round(e_gamma_0[i][0], 4),
        round(sd_gamma_0[i], 4)])
table.append(["h", round(e_h, 4), round(sd_h, 4), round(e_h_0, 4), round(sd_h_0, 4)])

print("\nApriorni a aposteriorni parametry:")
print(tabulate(table, headers, tablefmt="pretty"))

# ===== 6. Vypocet marginalni verohodnosti modelu metodou Gelfanda a Deye =====
# a) Zjednodusena varianta
gd_simple = 0
for i in range(s_1):
    gd_simple = gd_simple + 1/s_1 * 1/lik_ces(y, x, gamma[:,[i]], h[i])

print("\nMarginalni verohodnost CES produkcni funkce (zjednodusena varianta):")
print(f"Marginal likelihood = {round(1/gd_simple, 4)}")
print(f"Marginal likelihood (log) = {round(math.log(1/gd_simple), 4)}\n")

# b) plna varianta
theta = np.concatenate((gamma, np.array([h])))  # Slouceni matice gamma a vektoru h
theta_hat = t(np.array([np.mean(theta, axis=1)]))              # Sloupcovy vektor strednich hodnot
sigma_hat = np.cov(theta)
kk = theta_hat.size
pp = 0.01       # Pro (1 - p) procentni kvantil chi-kvadrat rozdeleni
chi_pp = chi2.ppf(1 - pp, kk)

# Simulace integracni konstanty pro omezenou apriorni hustotu
count_g = 0
for i in range(s):
    pom = gamma_0 + norm_rnd(v_0)   # p(gamma) ~ N(gamma_0,V_0)
    count_g += min(pom) > 0         # min(pom) > 0 vraci Bool - tedy True, nebo False, Python provadi konverzi na integer

int_c = 1/(count_g/s)[0]
gd = 0      # Statistika Gelfanda a Deye
for i in range(s_1):
    t_theta = t(theta[:,[i]] - theta_hat) @ inv(sigma_hat) @ (theta[:,[i]] - theta_hat)
    # Funkce hustoty vicerozmerneho rozdeleni
    f_theta = 1/(1 - pp) * 1/(2 * pi) ** (kk/2) * det(sigma_hat) ** (-1/2) * exp(1/2 * t_theta) if t_theta <= chi_pp else 0
    prior = prior_ces(gamma[:,[i]], h[i], gamma_0, v_0, h_0, nu_0)
    like = lik_ces(y, x, gamma[:,[i]], h[i])
    gd += f_theta/(int_c * prior * like) * 1/s_1
    progress_bar(i,s_1)

print("\nMarginalni verohodnost CES produkcni funkce (plna varianta):")
print(f"Marginal likelihood = {round(1/gd, 4)}")
print(f"Marginal likelihood (log) = {round(math.log(1/gd), 4)}\n")

# ===== 7. Srovnani skutecnych a modelovych momentu =====
e_y_ast = np.zeros((1, s_1))    # Vektor modelovych strednich hodnot
d_y_ast = np.zeros((1, s_1))    # Vektor modelovych rozptylu
std_y_ast = np.zeros((1, s_1))  # Vektor modelovych sm. odchylek

# Simulace (generovani umelych dat)
for i in range(s_1):
    f_xx = gamma[0, i] * (gamma[2, i] * x[:, [1]] ** gamma[3, i] + gamma[2, i] * x[:, [2]] ** gamma[3, i]) ** (1/gamma[3, i])
    y_ast = f_xx + uniform() * math.sqrt(1 / h[i])
    e_y_ast[0, i] = np.mean(y_ast)
    d_y_ast[0, i] = np.var(y_ast)
    std_y_ast[0, i] = np.std(y_ast)

mean_y = np.mean(y)
var_y = np.var(y)
std_y = np.std(y)
a = np.sum([val < mean_y for val in e_y_ast])
b = np.sum([val < var_y for val in d_y_ast])
c = np.sum([val < std_y for val in std_y_ast])

# Vypocet predikcnich (jednostrannych) p-hodnot
# a) pro stredni hodnotu
p_e = (a / s_1) * ((a / s_1) <= 0.5) + (1 - a / s_1) * ((a / s_1) > 0.5)
# b) pro rozptyl
p_d = (b / s_1) * ((b / s_1) <= 0.5) + (1 - b / s_1) * ((b / s_1) > 0.5)
# c) pro sm. odchylku
p_std = (c / s_1) * ((c / s_1) <= 0.5) + (1 - c / s_1) * ((c / s_1) > 0.5)

print("Predikcni p-hodnoty:")
print(f" p_E = {round(p_e, 4)}\n p_D = {round(p_d, 4)}\n p_std = {round(p_std, 4)}")

fig, axs = plt.subplots(3,1, tight_layout=True)

axs[0].hist(e_y_ast[0], bins=50, facecolor="blue", histtype="step")
axs[0].vline()
axs[1].hist(d_y_ast[0], bins=50, facecolor="blue", histtype="step")
axs[2].hist(std_y_ast[0], bins=50, facecolor="blue", histtype="step")

plt.show()