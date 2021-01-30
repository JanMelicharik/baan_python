# Funkce pro ziskani logaritmu jadra posteriorni hustoty

import numpy as np
from numpy import transpose as t
from numpy.linalg import inv

def log_post_ces(
    y: np.ndarray,
    x: np.ndarray,
    gamma: np.ndarray,
    h: float,
    gamma_0: np.ndarray,
    v_0: np.ndarray):

    if np.min(gamma) < 0:
        return -np.inf

    # CES funkce
    f_x = gamma[0][0] * (gamma[1][0] * x[:,[1]]**gamma[3][0] + gamma[2][0] * x[:,[2]]**gamma[3][0])**(1 / gamma[3][0])
    # Logaritmus (5.24) z Koop
    res = -1 / 2 * h * t(y - f_x) @ (y - f_x) -1 / 2 * t(gamma - gamma_0) @ inv(v_0) @ (gamma - gamma_0)

    return res[0][0]