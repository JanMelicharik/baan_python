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

data = pd.read_csv("energy_data.csv")
n = len(data.index)

# ===== 1. Priprava dat =====
# i + 1 protoze index jde v Pythonu od nuly
year = [i + 1 for i in data.index]

q = [val for val in data.log_quantity]
p = [val for val in data.log_price]
i = [val for val in data.income]
t = [val for val in data.date]

y = np.array([q]).T
x = np.array([[1]*len(q), year, p, [math.log(val) for val in i]]).T
n, k = x.shape

# ===== 2. Apriorni hustoty a apriorni hyperparametry =====
beta_0 = np.array(
    [
        [-30],
        [0.1],
        [-0.5],
        [0.5]
    ]
)
h_0 = 1 / 0.1 ** 2
nu_0 = 5
var_beta_0 = np.array(
    [
        [18 ** 2],
        [0.1 ** 2],
        [0.1 ** 2],
        [0.25 ** 2]
    ]
)
var_h_0 = 2 * h_0 ** 2 / nu_0
v_0 = np.diag(var_beta_0.ravel())   # Vstupem musi byt list nebo vektor

s_0 = 15_000 + 1
s_1 = 15_000
s = s_0 + s_1

beta = np.zeros((k, s))
h = np.zeros((1, s))
sd_nom = np.zeros((k, s))   # Citatel pro BF pomoci Savage-Dickey

# ===== 3. Gibbsuv vzorkovac =====
h[0] = h_0
beta[:, [0]] = beta_0
nu_1 = n + nu_0

for i in range(1, s):
    v_1 = inv(inv(v_0) + h[0][i - 1] * x.T @ x)
    beta_1 = v_1 @ (inv(v_0) @ beta_0 + h[0][i - 1] * x.T @ y)
    beta[:, [i]] = beta_1 + norm_rnd(v_1)
    h_1 = (1 / nu_1 * ((y - x @ beta[:,[i]]).T @ (y - x @ beta[:, [i]]) + nu_0 * h_0 ** (- 1))) ** (- 1)
    h[0][i] = gamm_rnd_koop2(h_1, nu_1)

    for j in range(k):
        # V Matlabu se jako skalovaci parametr vklada primo rozptyl, ktery je odmocnen ve funkci
        # V Pythonu se musi vlozit uz primo sm. odchylka
        sd_nom[j][i] = norm.pdf(0, beta_1[j][0], math.sqrt(v_1[j][j]))

    progress_bar(i, s)

# ===== 4. Vyhozeni prvnich s_0 vzorku =====
beta = beta[:, s_0:]
h = h[:, s_0:]
sd_nom = sd_nom[:, s_0:]

b_mean = np.mean(beta, axis=1)
h_mean = np.mean(h, axis=1)
b_var = np.var(beta, axis=1)
h_var = np.var(h, axis=1)

# ===== 5. Savage-Dickey pomer (Bayesuv faktor) =====
sd_denom = np.zeros((k, 1))
for i in range(k):
    sd_denom[i][0] = norm.pdf(0, beta_0[i][0], math.sqrt(v_0[i][i]))

bf = np.mean(sd_nom, axis=1) / sd_denom.ravel()

# ===== 6. Konvergencni diagnostiky (Geweke) =====
theta = np.c_[beta.T, h.T]
res_converg = geweke(theta)

# ===== 7. Prezentace vysledku =====
headers = ["", "Prior", "Posterior", "NSE", "Geweke CD", "Post. odds\n(beta_j=0)"]
table = []
for i in range(beta_0.shape[0]):
    table.append(
        [
            f"beta_{i+1}",
            f"{round(beta_0[i][0], 4)}\n({round(math.sqrt(var_beta_0[i][0]), 4)})",
            f"{round(b_mean[i], 4)}\n({round(math.sqrt(b_var[i]), 4)})",
            round(res_converg["nse"][i], 4),
            round(res_converg["cd"][i], 4),
            round(bf[i], 4)
        ]
    )

table.append(
    [
        "h",
        f"{round(h_0, 4)}\n({round(math.sqrt(var_h_0), 4)})",
        f"{round(h_mean[0], 4)}\n({round(math.sqrt(h_var[0]), 4)})",
        round(res_converg["nse"][4], 4),
        round(res_converg["cd"][4], 4),
        ""
    ]
)

print("\nApriorni a posteriorni parametry:")
print(tabulate(table, headers, tablefmt="pretty"))

fix, ax = plt.subplots(ncols=2, nrows=3, figsize=(10,10))
ax[0,0].set_title("beta_1")
ax[0,0].hist(beta[0, :], bins=30, histtype="bar", rwidth=0.9)
ax[0,1].set_title("beta_2")
ax[0,1].hist(beta[1, :], bins=30, histtype="bar", rwidth=0.9)
ax[1,0].set_title("beta_3")
ax[1,0].hist(beta[2, :], bins=30, histtype="bar", rwidth=0.9)
ax[1,1].set_title("beta_4")
ax[1,1].hist(beta[3, :], bins=30, histtype="bar", rwidth=0.9)
ax[2,0].set_title("h")
ax[2,0].hist(h[0, :], bins=30, histtype="bar", rwidth=0.9)
ax[2,1].remove()

plt.show()
# plt.show(block=False)
# pdb.set_trace()
