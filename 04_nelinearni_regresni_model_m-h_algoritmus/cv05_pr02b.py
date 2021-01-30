# Nastaveni cesty do domovske slozky
import sys
root_dir = "/".join(sys.argv[0].split("/")[:-3])
sys.path.insert(1, root_dir + "/Support/Python")

# importy pomocnych funkci
from progress_info import progress_bar

import numpy as np
import math
from tabulate import tabulate
import matplotlib.pyplot as plt

# ===== 1. Nastaveni independence chain M-H algoritmu =====

s = 50_000 + 1
s_0 = 30_000 + 1

theta = [0] * s
count = 0

# Kandidatska hustota N(theta[i-1], c**2)
c = 5

# ===== 2. Random walk chain M-H algoritmus =====

print("Random walk chain M-H algoritmus:")
for i in range(1, s):
    # a) Generovani kandidata theta_ast
    theta_ast = np.random.normal() * c + theta[i - 1]
    # b) Spocitani akceptacni pravdepodobnosti
    alpha = min(math.exp( -1/2 * abs(theta_ast) + 1/2 * abs(theta[i-1])), 1)
    # c) rozhodnuti o prijeti nebo odmitnuti kandidata
    if alpha > np.random.rand():
        theta[i] = theta_ast
        count += 1
    else:
        theta[i] = theta[i - 1]

    progress_bar(i, s)

# ===== 3. Prezentace vysledku =====
theta = theta[s_0:]     # Vyhozeni prvnich s_0 vzorku

# Graficke zobrazeni konvergence
fig, axs = plt.subplots(1,2, tight_layout=True)

axs[0].plot(theta[::500])
axs[0].set_title("Kazdy 500. odhad theta")
axs[1].hist(theta, bins=50, facecolor="blue", alpha=0.5, edgecolor="black")
axs[1].set_title("Rozdeleni theta")

# Pocitani momentu
e_theta = np.mean(theta)
d_theta = np.var(theta)
avg_count = count / (s - 1)     # Prumerna mira akceptace

headers = ["E[theta]", "D[theta]", "Prumerna akceptace"]
table = [[round(e_theta, 4), round(d_theta, 4), round(avg_count, 4)]]
print(tabulate(table, headers, tablefmt="pretty"))

plt.show()
