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

# Individual effects model - nehierarchicka apriorni hustota
# Apriorni hyperparametry
# Variabilni urovnove konstanty napric aerolinkami
alpha_0 = np.ones((1,n))
# Stejne ostatni parametry sklonu
beta_0 = np.concatenate((alpha_0, np.array([[1, 1, -1]])), axis=1).T
h_0 = 1 / 0.5 ** 2
var_alpha_0 = np.ones((1,n))
var_beta_0 = np.concatenate(
    (
        alpha_0,
        np.array([[1 ** 2, 1 ** 2, 1 ** 2]])
    ), axis=1).T
nu_0 = 1
var_h_0 = 2 * h_0 ** 2 / nu_0
v_0 = np.diag(var_beta_0)

mu_alpha_0 = 5
h_alpha_0 = 1 / (2 ** 2)
v_alpha_0 = 5
nu_alpha_0 = 1

s_0 = 1_000 + 1
s_1 = 10_000
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
