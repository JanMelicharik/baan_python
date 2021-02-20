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
y_i = np.reshape(y, (n,t))      # y_i pro i = 1 ... n (po sloupcich)

x = np.array([np.ones(n * t), data.Q, data.PF, data.LF]).T
x_i = np.reshape(x, (n,t,k))
x_i_til = np.reshape(x[:, 1:], (n,t,k-1))

x_ast = np.zeros((t * n, (n + k - 1)))
for i in range(n):
    x_ast[(t * i - t):(t * i), [i]] = np.ones((t, 1))
    x_ast[(t * i - t):(t * i), n:(n + k - 1)] = x_i_til[:, :, [i]]

# Random coefficients model - hierarchicka apriorni hustota
# Apriorni hyperparametry
beta_0 = np.zeros((k, n))
for i in range(n):
    beta_0[:, [i]] = np.array([[1, 1, 1, -1]]).T

h_0 = 1 / 0.5 ** 2
var_beta_0 = np.array([[1 ** 2, 1 ** 2, 1 ** 2, 1 ** 2]]).T
nu_0 = 1
var_h_0 = 2 * h_0 ** 2 / nu_0
v_0 = np.diag(var_beta_0)

mu_beta_0 = np.array([[1, 1, 1, -1]]).T
sigma_beta_0 = np.diag(var_beta_0)
v_beta_0 = v_0
nu_beta_0 = 1

s_0 = 1_000 + 1
s_1 = 10_000
s = s_0 + s_1

beta = np.zeros((k, s, n))
h = np.zeros((1, s))

v_beta = np.zeros((k, k, s))
mu_beta = np.zeros((k, s))

# Gibbsuv vzorkovac
h[:, [0]] = h_0
v_beta[:, :, [0]] = v_beta_0
mu_beta[:, [0]] = mu_beta_0

nu_beta_1 = n + nu_beta_0
nu_1 = t * n + nu_0

for i in range(s):
    for j in range(n):
        v_1 = inv(h[:, [i - 1]] * x_i[:, :, [j]].T @ x_i[:, :, [j]] + inv(v_beta[:, :, [i - 1]]))
        beta_1 = v_1 @ (
            h[:, [i - 1]] * x_i[:, :, [j]].T @ y_i[:, [j]] + \
            inv(v_beta[:, :, [i - 1]]) @ mu_beta[:, [i - 1]]
        )
        # Podminena hustota pro beta
        beta[:, [i], [j]] = beta_1 + norm_rnd(v_1)  # Koop 7.25
    
    # Podminena hustota pro mu_beta
    sigma_beta_1 = inv(n * inv(v_beta[:, :, [i - 1]]) + inv(sigma_beta_0))
    pom = 0
    for j in range(n):
        pom += beta[:, [i], [j]]

    mu_beta_1 = sigma_beta_1 @ (inv(v_beta[:, :, [i - 1]]) @ pom + inv(sigma_beta_0) @ mu_beta_0)
    mu_beta[:, [i]] = mu_beta_1 + norm_rnd(sigma_beta_1)

    pom = 0
    for j in range(n):
        pom += (beta[:, [i], [j]] - mu_beta[:, [i]] @ (beta[:, [i], [j]] - mu_beta[:, [i]]).T)

    v_beta_1 = pom + v_beta_0
    v_beta[:, :, [i]] = inv(wish_rnd(inv(nu_beta_1 @ v_beta_1), nu_beta_1))

    # Podminena hustota pro h
    pom = 0
    for j in range(n):
        pom += (
            y_i[:, [j]] - x_i[:, :, [j]] @ beta[:, [i], [j]]
        ).T @ (
            y_i[:, [j]] - x_i[:, :, [j]] @ beta[:, [i], [j]]
        )
    
    h_1 = (1 / nu_1 * (pom + nu_0 / h_0)) ** (- 1)
    h[:, [i]] = gamm_rnd_koop(h_1, nu_1)

# Vyhozeni prvnich s_0 vzorku
beta = beta[:, s_0:]
h = h[:, s_0:]
mu_beta = mu_beta[:, s_0:]

# Posteriorni stredni hodnoty
beta_mean = np.zeros((k, n))
for i in range(n):
    # Posteriorni stredni hodnota
    beta_mean[:, [i]] = np.mean(beta[:, :, i], axis=1)

h_mean = np.mean(h, axis=1)
mu_beta_mean = np.mean(mu_beta, axis=1)

beta_std = np.zeros((k, n))
for i in range(n):
    # Posteriorni sm. odchylka
    beta_std[:, [i]] = np.std(beta[:, :, [i]], axis=1)

h_std = np.std(h, axis=1)
mu_beta_std = np.std(mu_beta, axis=1)

#### PREZENTACE VYSLEDKU

