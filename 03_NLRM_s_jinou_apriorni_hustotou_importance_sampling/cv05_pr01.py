# Nastaveni cesty do domovske slozky
import sys
root_dir = "/".join(sys.argv[0].split("/")[:-3])
sys.path.insert(1, root_dir + "/Support/Python")

# importy pomocnych funkci
from progress_info import progress_bar
from my_boot import my_boot

# importy oficialnich knihoven
import numpy as np
from numpy import transpose as t
from math import exp
from tabulate import tabulate
import matplotlib.pyplot as plt

# ===== 1. Posteriorni hustota N(0,1) + Importance function U(a,b) =====
# Nastaveni importance sampler
s = 10_000
w = np.zeros((1,s))         # radkovy vektor vah
theta = np.zeros((1,s))     # vektor vzorku z importance function

# Importance function U(a,b)
a = -100
b = 100

# ===== 2. Importance sampling =====
for i in range(s):
    # a) Generovani vzorku z importance function, q(theta)
    theta[0,i] = np.random.rand() * (b - a) + a
    # b) Vypocet vah (na zaklade jadrovych hustot)
    p_theta = exp(-1 / 2 * theta[0,i]**2)
    q_theta = 1 / (b - a)
    w[0,i] = p_theta / q_theta
    progress_bar(i,s)

# ===== 3. Vypocet posteriorni stredni hodnoty a rozptylu =====

e_theta = w @ t(theta) / np.sum(w)          # vazena stredni hodnota
e_theta_2 = w @ (t(theta)**2) / np.sum(w)   # vazena stredni hodnota "na druhou"
d_theta = e_theta_2 - e_theta**2            # vazeny rozptyl
w_avg = np.mean(w)                          # prumerna vaha

headers = ["E[theta]", "D[theta]", "E[w]", "Vel. vzorku"]
table = [[round(e_theta[0][0], 4), round(d_theta[0][0], 4), round(w_avg, 4), s]]    # Musi byt list listu kvuli vytvoreni tabulky
print(f"\nImportance sampling s vyuzitim U({a},{b})")
print(tabulate(table, headers, tablefmt="pretty"))


# ===== 4. Vazeny bootstrap pro ziskani vyberu z posteriorni hodnoty =====
boot_smpl = my_boot(theta, s, w)
bootstrapped_theta = theta[0][boot_smpl]

fig, axs = plt.subplots(2,1, tight_layout=True)

axs[0].hist(theta, bins=50, facecolor="blue", alpha=0.5, edgecolor="black", histtype="step")
axs[1].hist(bootstrapped_theta, bins=50, facecolor="blue", histtype="step")

plt.show()