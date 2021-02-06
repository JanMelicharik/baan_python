# Funkce pro ziskani logaritmu jadra posteriorni hustoty

import numpy as np
from numpy.linalg import inv

def log_post_ces(y, x, gamma_ast, h, gamma_0, v_0):
    if np.min(gamma_ast) < 0:
        return -np.inf

    # CES funkce
    f_x = gamma_ast[0][0] * \
        ( gamma_ast[1][0] * x[:,[1]] ** gamma_ast[3][0] + gamma_ast[2][0] * x[:,[2]]**gamma_ast[3][0] ) ** \
        ( 1 / gamma_ast[3][0] )
    # Logaritmus (5.24) z Koop
    res = - h / 2 * (y - f_x).T @ (y - f_x) - 1 / 2 * (gamma_ast - gamma_0).T @ inv(v_0) @ (gamma_ast - gamma_0)
    return res
