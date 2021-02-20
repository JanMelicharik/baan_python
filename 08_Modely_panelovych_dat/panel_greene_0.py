import sys
sys.path.insert(1, "/".join(sys.path[0].split("/")[:-1]))

from support.progress_info import progress_bar
from support.norm_rnd import norm_rnd
from support.gamm_rnd_koop import gamm_rnd_koop

import pandas as pd
import numpy as np

import math
import pdb
import matplotlib.pyplot as plt
from tabulate import tabulate
from numpy.linalg import inv
from scipy.stats import norm

# Cost Data for U.S. Airlines, 90 Oservations On 6 Firms For 15 Years, 1970-1984
# Source: These data are a subset of a larger data set provided
# to the author by Professor Moshe Kim.
# They were originally constructed by Christensen Associates of Madison, Wisconsin.
#     * I = Airline,
#     * T = Year,
#     * Q = Output, in revenue passenger miles, index number,
#     * C = Total cost, in $1000,
#     * PF = Fuel price,
#     * LF = Load factor, the average capacity utilization of the fleet. 
data = pd.read_csv("data.csv")

n = 6   # Pocet aerolinek
t = 15  # Pocet let (1970 - 1984)
k = 4   # Pocet parametru

# Tvorba datovych matic
y = np.log(np.array([data.C]))  # Logaritmus nakladu
x = np.array([np.ones(n * t), data.Q, data.PF, data.LF]).T

# Pooled model
# "Standardni" Gibbsuv vzorkovac - jen vetsi matice a vektory
# Apriorni hyperparametry
beta_0 = np.array([[1, 1, 1, -1]]).T
h_0 = 1 / 0.5 ** 2
var_beta_0 = np.array([[1 ** 2, 1 ** 2, 1 ** 2, 1 ** 2]]).T 
nu_0 = 1
var_h_0 = 2 * h_0 ** 2 / nu_0
v_0 = np.diag(var_beta_0)

s_0 = 1_000 + 1
s_1 = 10_000
s = s_0 + s_1

beta = np.zeros((k, s))
h = np.zeros((1, s))

# Gibbsuv vzorkovac
h[0, [0]] = h_0
nu_1 = t * n + nu_0

for i in range(1, s):
    v_1 = inv(inv(v_0) + h[0, [i - 1]] * x.T @ x)
    beta_1 = v_1 @ (inv(v_0) @ beta_0 + h[0, [i - 1]] * x.T @ y)
    # Podminena hustota pro beta
    beta[:, [i]] = beta_1 + norm_rnd(v_1)
    h_1 = (1 / nu_1 * ((y - x @ beta[:, [i]]).T @ (y - x @ beta[:, [i]]) + nu_0 * 1 / h_0)) ** (- 1)
    # Podminena hustota pro h
    h[:, [i]] = gamm_rnd_koop(h_1, nu_1)

    progress_bar(i, s)

# Vyhozeni prvnich s_0 vzorku
beta = beta[:, s_0:]
h = h[:, s_0:]

beta_mean = np.mean(beta, axis=1)
h_mean = np.mean(beta, axis=1)
beta_std = np.std(beta, axis=1)
h_std = np.std(h, axis=1)

#### PREZENTACE VYSLEDKU



