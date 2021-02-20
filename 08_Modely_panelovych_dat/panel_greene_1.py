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

s_0 = 1_000 + 1
s_1 = 10_000
s = s_0 + s_1

beta = np.zeros((n + k - 1, s))
h = np.zeros((1, s))

# Gibbsuv vzorkovac
h[:, [0]] = h_0
nu_1 = t * n + nu_0

for i in range(1, s):
    v_1 = inv(inv(v_0) + h[0, [i - 1]] * x_ast.T @ x_ast)
    beta_1 = v_1 @ (inv(v_0) @ beta_0 + h[0, [i - 1]] @ x_ast.T @ y)
    # Podminena hustota pro beta
    beta[:, [i]] = beta_1 + norm_rnd(v_1)
    pom = 0
    for j in range(n):
        pom += (
            y_i[:, [j]] - beta[[j],[i]] * np.ones((t, 1)) - \
            x_i_til[:, :, [j]] @ beta[(n + 1):(n + k - 1), [i]]
        ).T @ (
            y_i[:, [j]] - beta[[j],[i]] * np.ones((t, 1)) - \
            x_i_til[:, :, [j]] @ beta[(n + 1):(n + k - 1), [i]]
        )

    h_1 = (1 / nu_1 * (pom + nu_0 / h_0)) ** (- 1)
    # Podminena hustota pro h
    h[:, [i]] = gamm_rnd_koop(h_1, nu_1)

# Vyhozeni prvnich s_0 vzorku
beta = beta[:, s_0:]
h = h[:, s_0:]

beta_mean = np.mean(beta, axis=1)
h_mean = np.mean(h, axis=1)
beta_std = np.std(beta, axis=1)
h_std = np.std(h, axis=1)

#### PREZENTACE VYSLEDKU
