import math
import numpy as np
from numpy.linalg import (inv, det)
from numpy.random import gamma as rnd_gamma

def prior_ces(gamma, h, gamma_0, v_0, h_0, nu_0):
    # Neomezena apriorni hustota (bez restrikci na kladnost gamma)
    k = gamma_0.shape[0]
    prior_gamma = 1 / (2 * math.pi) ** (k / 2) * det(v_0) ** (- 1 / 2) * \
                math.exp(- 1 / 2 * (gamma - gamma_0).T @ inv(v_0) @ (gamma - gamma_0))
    c_h = (2 * h_0 / nu_0) ** (nu_0 / 2) * math.gamma(nu_0 / 2)
    prior_h = 1 / c_h * h ** ((nu_0 - 2) / 2) * math.exp(- h * nu_0 / (2 * h_0))
    return prior_gamma * prior_h
