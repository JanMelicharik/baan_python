import sys
sys.path.insert(1, "/".join(sys.path[0].split("/")[:-1]))

from support.progress_info import progress_bar
from support.geweke import geweke
from support.norm_rnd import norm_rnd
from support.gamm_rnd_koop2 import gamm_rnd_koop2
from support.a_post import a_post

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
        [3.5],
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

z = np.array([year, p, [math.log(val) for val in i]]).T
p = z.shape[1]

beta = np.zeros((k, s))
h = np.zeros((1, s))
alpha_rw = np.zeros((p, s))

# po prvnich 50000 replikacich se vzal vysledny vektor rozptylu (a kovarianci) jako zaklad
# kovariancni matice; po te se dale ladi "d" pro ziskani zadouci
# akceptacni pravdepodobnosti;
# nebo - najde se modus posteriorni podminene hustoty pro alfa a prislusny
# hessian (podmineno apriornimi strednimi hodnotami, nebo jeste lepe
# posteriornimi strednimi hodnotami z prvnich 5000 replikaci

# postvar = np.eye(3)
# postvar = np.diag([0.132, 15.6423, 4.6199])
postvar = np.array(
    [
        [0.0775, -0.7646, 0.3380],
        [-0.7646, 16.5090, -8.191],
        [0.3380, -8.1914, 4.138]
    ]
)
d = 0.1
vscale_rw = d * postvar
alpha_0 = np.array(
    [
        [0.1278],
        [6.4524],
        [-3.5437]
    ]
)
# alpha_0 = np.array([[-0.1707], [-6.2956], [3.1309]])
# alpha_0 = np.array([[1.4918], [-7.4384], [4.3193]])

# ===== 3. Gibbsuv vzorkovac =====
h[0] = h_0
beta[:, [0]] = beta_0
alpha_rw[:, [0]] = alpha_0
nu_1 = n + nu_0
error_check = 0
count_rw = 0

for i in range(1, s):

    omega = np.diag(((z @ alpha_rw[:, i-1] + 1) ** 2).ravel())
    inv_omega = inv(omega)
    b_omega = inv(x.T @ inv_omega @ x) @ x.T @ inv_omega @ y
    v_1 = inv(inv(v_0) + h[0][i - 1] * x.T @ inv_omega @ x)
    beta_1 = v_1 @ (inv(v_0) @ beta_0 + h[0][i - 1] * x.T @ inv_omega @ x @ b_omega)
    try: 
        beta[:, [i]] = beta_1 + norm_rnd(v_1)
    except:
        beta[:, [i]] = beta_1 + norm_rnd(v_1 + 0.000000001)
        error_check += 1

    h_1 = (1 / nu_1 * \
            ((y - x @ beta[:, [i]]).T @ inv_omega @ (y - x @ beta[:, [i]]) + nu_0 * h_0 ** (- 1)) \
        ) ** (- 1)
    h[:, [i]] = gamm_rnd_koop2(h_1, nu_1)

    a_can_rw = alpha_rw[:, [i-1]] + norm_rnd(vscale_rw)
    log_accept_rw = min(
        a_post(
            a_can_rw,
            beta[:, [i]],
            h[0][i],
            y,
            x,
            z
        ) - a_post(
            alpha_rw[:, [i - 1]],
            beta[:, [i]],
            h[0][i],
            y,
            x,
            z
        ),
        0
    )

    if log_accept_rw > math.log(np.random.uniform()):
        alpha_rw[:, [i]] = a_can_rw
        count_rw += 1
    else:
        alpha_rw[:, [i]] = alpha_rw[:, [i - 1]]

    progress_bar(i, s)

# ===== 4. Vyhozeni prvnich s_0 vzorku =====
beta = beta[:, s_0:]
h = h[:, s_0:]
alpha_rw = alpha_rw[:, s_0:]

# Matice theta obsahuje parametry beta, h a alpha
theta = np.c_[beta.T, h.T, alpha_rw.T]
means = np.mean(theta, axis=0)
stds = np.std(theta, axis=0)

accept_ratio_rw = count_rw / s

# ===== 5. Konvergencni diagnostiky (Geweke) =====
res_converg = geweke(theta)

# HPDI pro beta a dalsi parametry
hpdi = np.quantile(theta, (0.05, 0.95), axis=0)

# ===== 7. Prezentace vysledku =====
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

fix, ax = plt.subplots(ncols=4, nrows=2, figsize=(12,8))
ax[0,0].set_title("beta_1")
ax[0,0].hist(theta[:, 0], bins=30, histtype="bar", rwidth=0.9)
ax[0,1].set_title("beta_2")
ax[0,1].hist(theta[:, 1], bins=30, histtype="bar", rwidth=0.9)
ax[0,2].set_title("beta_3")
ax[0,2].hist(theta[:, 2], bins=30, histtype="bar", rwidth=0.9)
ax[0,3].set_title("beta_4")
ax[0,3].hist(theta[:, 3], bins=30, histtype="bar", rwidth=0.9)
ax[1,0].set_title("alpha_1")
ax[1,0].hist(theta[:, 5], bins=30, histtype="bar", rwidth=0.9)
ax[1,1].set_title("alpha_2")
ax[1,1].hist(theta[:, 6], bins=30, histtype="bar", rwidth=0.9)
ax[1,2].set_title("alpha_3")
ax[1,2].hist(theta[:, 7], bins=30, histtype="bar", rwidth=0.9)
ax[1,3].set_title("h")
ax[1,3].hist(theta[:, 4], bins=30, histtype="bar", rwidth=0.9)


plt.show()
