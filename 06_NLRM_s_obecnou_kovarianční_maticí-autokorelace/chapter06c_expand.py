import sys
sys.path.insert(1, "/".join(sys.path[0].split("/")[:-1]))

from support.progress_info import progress_bar
from support.geweke import geweke
from support.norm_rnd import norm_rnd
from support.gamm_rnd_koop2 import gamm_rnd_koop2


import pandas as pd
import numpy as np

import math
import pdb
import matplotlib.pyplot as plt
from tabulate import tabulate
from numpy import transpose as t
from numpy.linalg import inv
from scipy.stats import norm

#######
# POZNAMKY:
# 1 - zmenit zobrazeni v tabulce - nazvy radku a sloupcu

# ===== 1. Priprava dat =====
# winning percentage in year t
# team on-base percentage in year t
# team slugging average in year t
# team earned run average in year t
data = pd.read_csv("yankees.tsv", sep="\t")

y_full = np.array([data.PCT]).T
x_full = np.array([np.ones(len(data.PCT)), data.OBP, data.SLG, data.ERA]).T

n, k = x_full.shape  # n - pocet pozorovani, k - pocet parametru

s_1 = 25_000
s_0 = 25_000 + 1
s = s_0 + s_1

# Neinformativni priory pro beta a h (nejsou uvedeny)
nu_0 = 0
inv_v_0 = np.zeros((k,k))

p = 3   # rad AR procesu
c = 1
rho_0 = np.zeros((p,1))
vrho_0 = c * np.eye(p)

beta = np.zeros((k, s))
h = np.zeros((1,s))
rho = np.zeros((p, s))

# Vytvoreni zpozdenych promennych + zkraceni puvodni rady
y = y_full[(p + 1):n, :]
x = x_full[(p + 1):n, :]

y_lag = np.zeros((n - p, p))
x_lag = np.zeros((n-p, k, p))
y_star = y
x_star = x

for i in range(p):
    y_lag[:, [i]] = y_full[(p + 1 - i):(n - i), [0]]
    x_lag[:, :, [i]] = x_full[(p + 1 - i):(n - i), :]
    y_star = y_star - rho_0[[i], 0] @ y_lag[:, [i]]
    x_star = x_star - rho_0[[i], 0] @ x_lag[:, :, [i]]

# ===== 2. Gibbsuv vzorkovac =====
# Pocatecni hodnoty
h[0][0] = 1
# pocet zamitnutych rho, lze pak vyuzit v ramci Savage-Dickey pro
# vypocet integracni konst. (pro p = 1 lze i analyticky)
count_post = 0
nu_1 = n
# Promenne S-D density ratio pro rho_i = 0
sd_prior = np.zeros((p, 1))
sd_post = np.zeros((p, s))

for i in range(1,s):
    y_star = y
    x_star = x

    for j in range(p):
        y_star = y_star - rho[j][s - 1] * y_lag[:, j]
        x_star = x_star - rho[j][s - 1] * x_lag[:, :, [j]] 

    v_1 = inv(h[0][i - 1] * x_star.T @ x_star)
    beta_1 = v_1 @ (h[0][i - 1] * x_star.T @ x_star.T)
    beta[:, [i]] = beta_1 + norm_rnd(v_1)

    h_1 = (
        1 / nu_1 * (
            (y_star - x_star @ beta[:, [i]]).T @ (y_star - x_star @ beta[:, [i]])
        )
    ) ** (- 1)
    h[0][i] = gamm_rnd_koop2(h_1, nu_1)

    # Vypocet a sestrojeni matice E (nahodnych slozek a jejich zpozdeni)
    eps_full = y_full - x_full @ beta[:, [i]]
    eps = eps_full[p:n, :]
    e = np.zeros((n - p, p))
    for j in range(p):
        e[:, [j]] = eps_full[(p + 1 - j):(n - j), [0]]

    vrho_1 = inv(inv(vrho_0) + h[0][i] * e.T @ e)
    rho_1 = vrho_1 @ (inv(vrho_0) @ rho_0 + h[0][i] * e.T @ eps)

    # Vyber rho jen ze stacionarni oblasti
    test_rho = 0
    while test_rho == 0:
        rho[:, [i]] = rho_1 + norm_rnd(vrho_1)
        count_post += 1
        # ZDE JE POTREBA PRIDAT flipud() alternativu

        if abs(rho[:, [i]]) < 1:    # MOJE POZNAMKA: all() mozna?
            test_rho = 1
            count_post -= 1

    for j in range(p):
        sd_post[j][i] = norm.pdf(0, rho[j][i], vrho_1[j][j])

    progress_bar(i, s)

# Vyhozeni prvnich s_0 vzorku
beta = beta[:, s_0:]
h = h[:, s_0:]
rho = rho[:, s_0:]
sd_post = sd_post[:, s_0:]

# Vypocet posteriornich charakteristik
theta = np.c_[beta.T, h.T, rho.T]
means = np.mean(theta, axis=0)
stds = np.std(theta, axis=0)

# ===== 3. Konvergencni diagnostiky =====
res_converg = geweke(theta)
hpdi = np.quantile(ksi, (0.05, 0.95), axis=0)

# ===== 4. Savage-Dickey =====
# intergracni konst. post. hustoty
sd_post_int = 1 - (count_post / (s + count_post))
sd_nom = np.mean(sd_post, axis=1)
# Savage-Dickey prior
for i in range(p):
    sd_prior[i][0] = norm.pdf(0, rho_0[i][0], vrho_0[i][i])

# integracni konstanta pro apriorni hustotu
# pocet zamitnutych chyb
count_prior = 0
for i in range(s):
    # vyber rho pouze ze stacionarni oblasti
    test_rho = True
    while test_rho:
        rho_prior = rho_0 + norm_rnd(vrho_0)
        count_prior += 1
        # ZDE JE POTREBA PRIDAT flipud() alternativu

# Integracni konst. apr. hustoty
sd_prior_int = 1 - (count_prior / (s + count_prior))
# Bayesuv faktor rho_i = 0
bf = sd_nom / sd_prior @ (sd_prior_int / sd_post_int)

# ===== 5. Prezentace vysledku =====
headers = ["", "Prior", "Posterior", "NSE", "Geweke CD", "95% HPDI"]
params = ["beta_1","beta_2","beta_3","beta_4","h","alpha_1","alpha_2","alpha_3"]
table = []
for i in range(theta.shape[1]):
    table.append(
        [
            params[i],
            "",
            f"{round(means[i], 4)}\n({round(stds[i], 4)})",
            round(res_converg["nse"][i], 4),
            round(res_converg["cd"][i], 4),
            f"[{round(hpdi[0][i], 4)}; {round(hpdi[1][i], 4)}]"
        ]
    )

for i in range(beta_0.shape[0]):
    table[i][1] = f"{round(beta_0[i][0], 4)}\n({round(math.sqrt(var_beta_0[i][0]), 4)})"

table[4][1] = f"{round(h_0, 4)}\n({round(math.sqrt(var_h_0),4)})"
for i in range(5,8):
    table[i][1] = f"neinf.\n(N.A.)"

print("\nApriorni a posteriorni parametry:")
print(tabulate(table, headers, tablefmt="pretty"))
print(
    f"Pocet replikaci:"\
    f"\n - {s} celkem"\
    f"\n - {s_0} vyrazeno"\
    f"\nMira akceptace: {round(accept_ratio_rw, 4)}"\
    f"\nChyb: {error_check}"
)
