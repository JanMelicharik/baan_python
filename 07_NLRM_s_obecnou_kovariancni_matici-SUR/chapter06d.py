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
data_yankees = pd.read_csv("yankees.tsv", sep="\t")
data_redsox = pd.read_csv("redsox.tsv", sep="\t")

# SUR model - tvorba struktury
m = 2
n = len(data_yankees.PCT)

y_y = np.array([data_yankees.PCT]).T
r_y = np.array([data_redsox.PCT]).T

y_x = np.array(
    [
        np.ones(len(data_yankees.PCT)),
        data_yankees.OBP,
        data_yankees.SLG,
        data_yankees.ERA
    ]
).T

r_x = np.array(
    [
        np.ones(len(data_redsox.PCT)),
        data_redsox.OBP,
        data_redsox.SLG,
        data_redsox.ERA
    ]
).T

k_1 = y_x.shape[1]
k_2 = r_x.shape[1]
k = k_1 + k_2

x = np.zeros((m * n, k))
y = np.zeros((m * n, 1))
for i in range(n):
    x[(m * i - 1), 0:k_1] = y_x[i, :]
    x[(m * i), k_1:k] = r_x[i, :]
    y[(m * i - 1), :] = y_y[i,:]
    y[(m * i), :] = r_y[i,:]

# Apriorni hyperparametry a nastaveni Gibbsova vzorkovace
s_0 = 15_000 + 1    # Pocatecni podminka
s_1 = 15_000
s = s_0 + s_1

# Neinformativni priory pro beta a h (nejsou uvedeny)
beta_0 = np.zeros((k, 1))
v_0 = 4 * np.eye(k)
nu_0 = 0
inv_h_0 = np.zeros((m,m))
h_0 = np.eye(m)

# trojrozmerna matice presnosti chyb (treti rozmer = beh Gibbsova vzorkovace)
beta = np.zeros((k, s))
h = np.zeros((1,s))
# Koop 6.52
nu_1 = nu_0 + n

beta[:, [0]] = beta_0
h[:, :, [0]] = h_0

# ===== 2. Gibbsuv vzorkovac =====
for i in range(1,s):
    sum_v_1 = np.zeros((k, k))
    sum_b_1 = np.zeros((k, k))
    for j in range(n):
        sum_v_1 += x[(j * 2):(j * 2 + 2), :].T @ h[:, :, [j - 1]] @ x[(j * 2):(j * 2 + 2), :]
        sum_b_1 += x[(j * 2):(j * 2 + 2), :].T @ h[:, :, [j - 1]] @ y[(j * 2):(j * 2 + 2), [0]]

    v_1 = inv(v_0) + sum_v_1    # Koop 6.50
    beta_1 = inv(v_1) @ (inv(v_0) @ beta_0 + sum_b_1)   # Koop 6.51

    # Podminena hustota pro beta
    beta[:, [i]] = beta_1 + norm_rnd(inv(v_1))  # Koop 6.49
    # Podminena hustota pro h
    sum_h = np.zeros((m, m))
    for j in range(n):
        # cast sumy v Koop 6.54
        sum_h += (
                y[(j * 2):(j * 2 + 2), [0]] - \
                x[(j * 2):(j * 2 + 2), :] @ beta[:, [i]]
            ) @ (
                y[(j * 2):(j * 2 + 2), [0]] - \
                x[(j * 2):(j * 2 + 2), :] @ beta[:, [i]]
            ).T

    # Koop 6.52
    h_1 = inv(inv_h_0 + sum_h)
    # Podminena hustota pro H
    # Koop 6.49
    h[:, :, [i]] = wish_rnd(h_1, nu_1)

    progress_bar(i, s)

# Vyhozeni prvnich s_0 vzorku
beta = beta[:, s_0:]
h = h[:, s_0:]
# Korelace nahodne slozky napric rovnicemi
r = h[1, 0, :] / (math.sqrt(h[0, 0, :] * h[1, 1, :]))

# Vypocet posteriornich charakteristik
beta_mean = np.mean(beta, axis=1)
h_mean = np.mean(h, axis=1)
r_mean = np.mean(r, axis=1)
beta_std = np.std(beta, axis=1)
h_std = np.std(r, axis=1)
r_std = np.std(h, axis=1)

hpdi_beta = np.quantile(beta, (0.05, 0.95), axis=1)
hpdi_r = np.quantile(r, (0.05, 0.95), axis=1)

##### ZBYVA DODELAT


# ===== 4. Prezentace vysledku =====
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
