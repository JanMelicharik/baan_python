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

# country     country code
# year        olympics year
# gdp         gross domestic product, 1995 dollars
# pop         population
# gold        number of gold medals won
# silver      number of silver medals won
# bronze      number of bronze medals won
# medaltot    total number of medals won
# host        = 1 if host country
# planned     = 1 if non-soviet planned
# soviet      = 1 if soviet
data_raw = pd.read_csv("olympics.csv", sep=";", decimal=",")
data = pd.DataFrame(data_raw.share, data_raw.gdp, data_raw.pop).dropna()

# Priprava dat - prace s datovym souborem jiz ocistenym o NA
y = np.array([data.share]).T
x = np.array(
    [
        np.ones((1, len(y))),
        np.log(data.gdp),
        np.log(data.pop)
    ]
).T

# Nastaveni apriornich hyperparametru a Gibbsova vzorkovace
# Apriorni hyperparametry
# p(beta) ~ N(beta_0, v_0)
# p(h) ~ Gamma(h_0, nu_0)
beta_0 = np.array([[0, 0, 0]]).T
v_0 = np.diag([0.2 ** 2, 0.01 ** 2, 0.01 ** 2])
nu_0 = 100
s2_0 = 0.1 ** 2
h_0 = 1 / s2_0

# Nastaveni Gibbsova vzorkovace
s_0 = 30_000 + 1
s_1 = 20_000
s = s_0 + s_1

beta = np.zeros((len(beta_0), s))
h = np.zeros((1, s))
y_ast = np.zeros((len(y), s))

n = len(y)
nu_1 = nu_0 + n
# Gibbsuv vzorkovac
for i in range(1, s):
    # 1. Blok - podm. hustota p(beta|h, y_ast) ~ N(beta_1, v_1)
    # Koop 4.4
    v_1 = inv(inv(v_0) + h[0, [i - 1]] * (x.T @ x))
    # Koop 4.5
    beta_1 = v_1 @ (inv(v_0) @ beta_0 + h[0, [i - 1]] * (x.T @ y_ast[:, [i - 1]]))
    # Koop 4.7
    beta[:, [i]] = beta_1 + norm_rnd(v_1)

    # 2. Blok - podm. hustota p(h|beta, y_ast) ~ Gamma(h_1, nu_1)
    # Koop 4.10
    h_1 = (1 / nu_1 * \
        (y_ast[:, [i - 1]] - x @ beta[:, [i]]).T @ (y_ast[:, [i - 1]] - x @ beta[:, [i]]) + \
    nu_0 / h_0) ** (- 1)
    # Koop 4.8
    h[:, [i]] = gamm_rnd_koop(h_1, nu_1)

    # 3. Blok - podm. hustota pro p(y_ast|beta, h, y)
    y_ast[:, [i]] = y
    

