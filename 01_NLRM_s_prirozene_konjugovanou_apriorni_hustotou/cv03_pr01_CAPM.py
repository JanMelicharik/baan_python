# Nastaveni cesty do domovske slozky
import sys
root_dir = "/".join(sys.argv[0].split("/")[:-2])
sys.path.insert(1, root_dir + "/support")

# Importy podpurnych funkci ze slozky Support
from my_nlrm import my_nlrm
from gamm_rnd_koop import gamm_rnd_koop
from norm_rnd import norm_rnd
from progress_info import progress_bar

from tabulate import tabulate
from math import exp
from scipy.stats import t as student
from numpy import ones,\
                  zeros,\
                  transpose as t

import pandas as pd
import numpy as np
import math

# ! Poznamka ke scriptu: !
# promenna beta_0 a cov_beta_0 musi byt zadane jako prvky 2D matice
# i v pripade, ze jde pouze o skalary (cili, vektory delky 1).
# Prvky 2D matice se zapisuji jako list listu --> [[1]]

# Dve moznosti nacteni dat (u CSV souboru je nutne specifikovat
# oddelovac sloupcu a nastavit znak pro desetinnou carku).
data = pd.read_excel(root_dir + "/data/capm2_data.xlsx")
# data = pd.read_csv(root_dir + "/data/capm2_data.csv", delimiter=";", decimal=",")

# ===== 1. Odhad modelu =====
y_0 = pd.DataFrame({"y": data.GM - data.RKFREE})
x_0 = pd.DataFrame({"const": 1, "x": data.MKT - data.RKFREE})

y_0 = np.array(y_0)
x_0 = np.array(x_0)

# Nastaveni apriornich hyperparametru
# beta|h ~ N(beta_0,V_0)
# h ~ G(h_0,nu_0)

s2_0 = 0.2 ** 2     # skalar
h_0 = 1 / s2_0      # skalar
nu_0 = 10           # skalar

beta_0 = np.array([[0],[1]])                    # vektor (matice 2x1)
cov_beta_0 = np.diag([0.05 ** 2, 0.5 ** 2])     # matice
v_0 = cov_beta_0 * (nu_0 - 2) / nu_0 * h_0      # matice

res_gm = my_nlrm(y_0, x_0, beta_0, v_0, h_0, nu_0)

headers = ["Parametr", "Prior", "Prior std.", "Posterior", "Posterior std."]
table_1 = [
    [
        "Alpha", 
        res_gm["beta_0"][0][0],             # prvni [0] bere prvni prvek vektoru
        res_gm["b0_std"][0][0],             # druha [0] rozbaluje list
        round(res_gm["beta_1"][0][0], 4),   # (jde pouze o vizualni stranku)
        round(res_gm["b1_std"][0][0], 4)    # (odstranuje hranate zavorky)
    ],
    [
        "Beta",
        res_gm["beta_0"][1][0],
        res_gm["b0_std"][1][0],
        round(res_gm["beta_1"][1][0], 4),
        round(res_gm["b1_std"][1][0], 4)
    ],
    [
        "h",
        round(res_gm["h_0"], 4),
        round(res_gm["h0_std"], 4),
        round(res_gm["h_1"][0][0], 4),
        round(res_gm["h1_std"], 4)
    ]
]

print("Odhad NLRM s NCP:")
print(tabulate(table_1, headers, tablefmt="pretty"), "\n")

# ===== 2. Test hypotezy, ze beta = 1 =====
y_1 = np.array(pd.DataFrame({"y": data.GM - data.MKT}))
x_1 = ones(y_1.shape)

# vstupy do funkce my_nlrm musi byt matice, viz vyse
beta_0_rest = np.array([beta_0[0]])
v_0_rest = np.array([[v_0[0][0]]])

res_gm_rest = my_nlrm(y_1, x_1, beta_0_rest, v_0_rest, h_0, nu_0)

table_2 = [
    [
        "Alpha", 
        res_gm_rest["beta_0"][0][0],
        res_gm_rest["b0_std"][0][0],
        round(res_gm_rest["beta_1"][0][0], 4),
        round(res_gm_rest["b1_std"][0][0], 4)
    ],
    [
        "h",
        round(res_gm_rest["h_0"], 4),
        round(res_gm_rest["h0_std"], 4),
        round(res_gm_rest["h_1"][0][0], 4),
        round(res_gm_rest["h1_std"], 4)
    ]
]

print("Omezeny model (beta = 1):")
print(tabulate(table_2, headers, tablefmt="pretty"), "\n")

log_bf = res_gm_rest["log_ml"] - res_gm["log_ml"]
bf = exp(log_bf)    # Bayesuv faktor (odlogaritmujeme predchozi vyraz)

print("Bayesuv faktor porovnavajici omezeny a neomezeny model:")
print(f"BF = {round(bf, 4)}", "\n")

# ===== 3. Hypoteza, ze beta > 1 =====
mc = 100_000 # pocet simulaci      (zvyste napr. na 10_000)
beta_sim = np.array([[],[]])

print("Vypocet pravd., ze beta > 1 pomoci simulace:")
for i in range(mc):
    h_sim = float(gamm_rnd_koop(res_gm["h_1"], res_gm["nu_1"], (1,1)))
    new_column = norm_rnd(1/h_sim * res_gm["v_1"]) + res_gm["beta_1"]
    beta_sim = np.append(beta_sim, new_column, axis=1)
    progress_bar(i, mc)

# Vypocet pravd. beta > 1
pr_beta = sum(t(beta_sim > 1))[1] / mc
print(f"Pravdepodobnost, ze beta > 1:")
print(f"Pr. = {round(pr_beta, 4)}")

# Analyticky vypocet pravdepodobnost
# a) standardizace skalovaneho t-rozdeleni (p(beta|y)) pro druhy prvek vektoru parametru beta
zscore = float((1 - res_gm["beta_1"][1]) / res_gm["b1_std"][1])

# b) vypocet odpovidajiciho kvantilu ze standardizovaneho centrovaneho t-rozdeleni
pr_beta_analyticky = 1 - student.cdf(zscore, res_gm["nu_1"])
print(f"Pr. = {round(pr_beta_analyticky, 4)} (analyticky)")
