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

beta = np.zeros((k, s))
h = np.zeros((1,s))

# ===== 2. Gibbsuv vzorkovac =====
# Pocatecni hodnoty
h[0][0] = 1
# pocet zamitnutych rho, lze pak vyuzit v ramci Savage-Dickey pro
# vypocet integracni konst. (pro p = 1 lze i analyticky)
count_post = 0
nu_1 = n

for i in range(1,s):
    v_1 = inv(h[0][i-1] * x.T @ x)
    beta_1 = v_1 @ (h[0][i - 1] * x.T @ y)
    beta[:, [i]] = beta_1 + norm_rnd(v_1)

    h_1 = (1 / nu_1 * ((y - x @ beta[:, [i]]).T @ (y - x @ beta[:, [i]]))) ** (- 1)
    h[0][i] = gamm_rnd_koop2(h_1, nu_1)

    progress_bar(i, s)

# Vyhozeni prvnich s_0 vzorku
beta = beta[:, s_0:]
h = h[:, s_0:]

# Vypocet posteriornich charakteristik
theta = np.c_[beta.T, h.T]
means = np.mean(theta, axis=0)
stds = np.std(theta, axis=0)

# ===== 3. Konvergencni diagnostiky =====
res_converg = geweke(theta)
hpdi = np.quantile(ksi, (0.05, 0.95), axis=0)

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
