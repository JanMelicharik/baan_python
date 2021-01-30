# Verohodnostni funkcep pro CES produkcni funkci

import numpy as np
from numpy import transpose as t
import math

def lik_ces(
    y: np.ndarray,
    x: np.ndarray,
    gamma: np.ndarray,
    h: float):

    f_x = gamma[0][0] * (gamma[1][0] * x[:,[1]]**gamma[3][0] + gamma[2][0] * x[:,[2]]**gamma[3][0])**(1 / gamma[3][0])
    n = y.shape[0]
    log_res = -n/2 * math.log(2 * math.pi) + n/2 * math.log(h) - h/2 * t(y - f_x) @ (y - f_x)
    res = math.exp(log_res)
    
    return res