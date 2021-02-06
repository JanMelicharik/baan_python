# Verohodnostni funkce pro CES produkcni funkci
import numpy as np
import math

def lik_ces(y, x, gamma, h):

    f_x = gamma[0][0] * \
        ( gamma[1][0] * x[:,[1]] ** gamma[3][0] + gamma[2][0] * x[:,[2]] ** gamma[3][0] ) ** \
        ( 1 / gamma[3][0] )
    n = y.shape[0]
    log_res = - n / 2 * math.log(2 * math.pi) + \
                n / 2 * math.log(h) - h / 2 * \
                (y - f_x).T @ (y - f_x)
    return math.exp(log_res)
