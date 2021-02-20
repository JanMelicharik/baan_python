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

# Source: Spector and Mazzeo (1980).
#   * Obs = observation,
#   * TUCE = Test score on economics test,
#   * PSI = participation in program,
#   * GRADE = Grade increase (1) or decrease (0) indicator 
#   * GPA = grade point average
data = pd.read_csv("data.csv")

k = 4   # Pocet parametru
n = 32  # Pocet pozorovani

y = np.array([data.GRADE])
x = np.array([np.ones(len(data.GRADE)), data.GPA, data.TUCE, data.PSI])

# Apriorni hyperparametry
beta_0 = np.array([[0, 0, 0, 0]]).T
nu_0 = 2

# Pro identifikovatelnost parametru volime primo h = 1
h = 1

var_beta_0 = np.array([[10 ** 2, 10 ** 2, 10 ** 2, 10 ** 2]]).T
v_0 = np.diag(var_beta_0)

y_ast_0 = y

s_0 = 5_000 + 1
s_1 = 10_000
s = s_0 + s_1

beta = np.zeros((k, s))
y_ast = np.zeros((n, s))


