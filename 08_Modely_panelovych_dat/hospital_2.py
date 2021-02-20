import sys
sys.path.insert(1, "/".join(sys.path[0].split("/")[:-1]))

# importy pomocnych funkci
from norm_rnd import norm_rnd
from gamm_rnd_koop import gamm_rnd_koop
from progress_info import progress_bar
from geweke import geweke

# importy oficialnich knihoven
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tabulate import tabulate
from scipy.io import loadmat
from numpy import transpose as t
from numpy.linalg import inv
from math import sqrt
from numpy import diag,\
                  ones,\
                  zeros

raw_data = loadmat("data_hospital.mat")

# Data jsou dictionary, kde hodnota je list hodnot
variables = [
    "hospital_id",  # ... hospital ID
    "year",         # ... year
    "costs",        # ... log of hospital operating costs (in thousands of dollars)
    "beds",         # ... log of number of beds in hospital
    "inpatient",    # ... log of number of inpatients visit
    "outpatient",   # ... log of number of outpatient visit
    "case_mix",     # ... log of case mix index (higher values = difficult cases)
    "k",            # ... log of capital stock (in thousand dollars)
    "nonprofit",    # ... 1 for non-profit hospitals (=0 otherwise)
    "forprofit",    # ... 1 for for-profit hospitals (=0 otherwise)
]
data = {}
for variable in variables:
    data[variable] = raw_data[variable].transpose().tolist()[0]

n = 382     # pocet nemocnic
t = 5       # pocet let
k = 8       # pocet parametru

# Tvorba datovych matic
y = np.array([data["costs"]])
y_i = np.reshape(y, (n,t))

x = np.array([
        np.ones(n * t),
        data["beds"],
        data["inpatient"],
        data["outpatient"],
        data["case_mix"],
        data["k"],
        data["nonprofit"],
        data["forprofit"]
    ])
# Zkontrolovat, ze data jsou spravne usporadana
x_i = np.reshape(x, (n,t,k))
x_i_til = np.reshape(x[:, 1:], (n,t,k-1))

x_ast = np.zeros((t * n, (n + k - 1)))
for i in range(n):
    x_ast[(t * i - t):(t * i), [i]] = np.ones((t, 1))
    x_ast[(t * i - t):(t * i), n:(n + k - 1)] = x_i_til[:, :, [i]]

# Individual effects model - nehierarchicka apriorni hustota
# Apriorni hyperparametry
# Variabilni urovnove konstanty napric aerolinkami
alpha_0 = np.ones((1,n))
# Stejne ostatni parametry sklonu
beta_0 = np.array([[1, 1, 1, 1, 1, 1, 1]])
h_0 = 1 / 1 ** 2
var_alpha_0 = np.ones((1,n))
var_beta_0 = np.array([[1 ** 2, 1 ** 2, 1 ** 2, 1 ** 2, 1 ** 2, 1 ** 2, 1 ** 2]])
nu_0 = 1
var_h_0 = 2 * h_0 ** 2 / nu_0
v_0 = np.diag(var_beta_0)

mu_alpha_0 = 5
h_alpha_0 = 1 / (2 ** 2)
v_alpha_0 = 5
nu_alpha_0 = 1

s_0 = 500 + 1
s_1 = 500
s = s_0 + s_1

beta = np.zeros((k - 1, s))
h = np.zeros((1, s))
alpha = np.zeros((n, s))

v_alpha = np.zeros((1, s))
mu_alpha = np.zeros((1, s))

# Gibbsuv vzorkovac
h[:, [0]] = h_0
v_alpha[:, [0]] = v_alpha_0
mu_alpha[:, [0]] = mu_alpha_0
alpha[:, [0]] = alpha_0.T

nu_1  = t * n + nu_0
nu_alpha_1 = nu_alpha_0 + n

for i in range(1, s):
    pom = 0
    for j in range(n):
        pom += x_i_til[:, :, [j]].T @ x_i_til[:, :, [j]]

    v_1 = inv(inv(v_0) + h[:, [i - 1]] * pom)
    
    pom = 0
    for j in range(n):
        pom += x_i_til[:, :, [j]].T @ (y_i[:, [j]] - alpha[j, [i - 1]] @ np.ones((t, 1)))
    # Podminena hustota pro beta
    beta_1 = v_1 @ (inv(v_0) @ beta_0 + h[:, [i - 1]] * pom)

    pom = 0
    for j in range(n):
        pom += (
            y_i[:, [j]] - alpha[j, [i - 1]] @ np.ones((t, 1)) - x_i_til[:, :, [j]] @ beta[:, i]
        ).T @ (
            y_i[:, [j]] - alpha[j, [i - 1]] @ np.ones((t, 1)) - x_i_til[:, :, [j]] @ beta[:, i] 
        )
    # Podminena hustota pro h
    h_1 = (1 / nu_1 * (pom + nu_0 / h_0)) ** (- 1)
    h[:, [i]] = gamm_rnd_koop(h_1, nu_1)

    # Podminene hustoty pro kazde alpha_j
    for j in range(n):
        v_1_j = (v_alpha[:, [i - 1]] * h[:, [i]] ** (- 1)) / (t * v_alpha[:, [i - 1]] + h[:, [i]] ** (- 1))
        alpha_1 = (
            v_alpha[:, [i - 1]] @ (y_i[:, [j]] - x_i_til[:, :, [j]] @ beta[:, [j]]).T @ \
            np.ones((t, 1)) + h[:, [i]] ** (- 1) * mu_alpha[:, [i - 1]]
        ) / (
            t * v_alpha[:, [i - 1]] + h[:, [i]] ** (- 1)
        )
        alpha[[j],[s]] = alpha_1 + norm_rnd(v_1_j)      # Koop 7.16

    # Podminena hustota pro mu_alpha
    sigma_2_alpha_1 = v_alpha[:, [i - 1]] / h_alpha_0 / (v_alpha[:, [i - 1]] + n / h_alpha_0)
    mu_alpha_1 = (
        v_alpha[:, [i - 1]] * mu_alpha_0 / h_alpha_0 * np.sum(alpha[:, [i]])
    ) / (
        v_alpha[:, [i - 1]] + n / h_alpha_0
    ) 
    mu_alpha[:, [i]] = mu_alpha_1 + norm_rnd(sigma_2_alpha_1)   # Koop 7.17

    # Podminena hustota pro inv(v_alpha)
    pom = 0
    for j in range(n):
        pom += (alpha[[j], [i]] - mu_alpha[:, [i]]) ** 2

    v_alpha_1 = (pom + v_alpha_0 * nu_alpha_0) / nu_alpha_1
    v_alpha[:, [i]] = inv(gamm_rnd_koop(inv(v_alpha_1), nu_alpha_1))    # Koop 7.18

# Vyhozeni prvnich s_0 vzorku
alpha = alpha[:, s_0:]
beta = beta[:, s_0:]
h = h[:, s_0:]

# Posteriorni stredni hodnoty
alpha_mean = np.mean(alpha, axis=1)
beta_mean = np.mean(beta, axis=1)
h_mean = np.mean(h, axis=1)
# Posteriorni sm. odchylky
alpha_std = np.std(alpha, axis=1)
beta_std = np.std(beta, axis=1)
h_std = np.std(h, axis=1)

#### PREZENTACE VYSLEDKU
