# Vraci inverzni hodnotu cdf Gamma rozdeleni pro danou hodnotu p
# Pro uplny popis viz. tuto funkci pro MATLAB

import math

x = max(a - 1, 0.1)
dx = 1
eps = 2**(-52)

while abs(dx) > 256 * eps * max(x, 1)

