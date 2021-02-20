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

# Pooled model
# "standardni" Gibbsuv vzorkovac - jen vetsi matice a vektory
# Apriorni hyperparametry

beta_0 = np.array([[1, 1, 1, 1, 1, 1, 1, 1]]).T
h_0 = 1 / 0.5 ** 2
var_beta_0 = np.array([[1 ** 2, 1 ** 2, 1 ** 2, 1 ** 2, 1 ** 2, 1 ** 2, 1 ** 2, 1 ** 2]]).T
nu_0 = 1
var_h_0 = 2 * h_0 ** 2 / nu_0
v_0 = np.diag(var_beta_0)

s_0 = 1_000 + 1
s_1 = 10_000
s = s_0 + s_1

beta = np.zeros((k, s))
h = np.zeros((1, s))

# Gibbsuv vzorkovac
h[:, [0]] = h_0
nu_1 = t * n + nu_0

for i in range(1, s):
    v_1 = inv(inv(v_0) + h[0, [i - 1]] * x.T @ x)
    beta_1 = v_1 @ (inv(v_0) * beta_0 + h[0, [i - 1]] * x.T @ y)
    # Podminena hustota pro beta
    beta[:, [i]] = beta_1 + norm_rnd(v_1)
    h_1 = (1 / nu_1 * ((y - x @ beta[:, [i]]).T @ (y - x @ beta[:, [i]]) + nu_0 / h_0)) ** (- 1)
    # Podminena hustota pro h
    h[:, [i]] = gamm_rnd_koop(h_1, nu_1)

# Vyhozeni prvnich s_0 vzorku
beta = beta[:, s_0:]
h = h[:, s_0:]

beta_mean = np.mean(beta, axis=1)
h_mean = np.mean(h, axis=1)
beta_std = np.std(beta, axis=1)
h_std = np.std(h, axis=1)

# Konvergencni diagnostiky se neresi

### PREZENTACE VYSLEDKU
