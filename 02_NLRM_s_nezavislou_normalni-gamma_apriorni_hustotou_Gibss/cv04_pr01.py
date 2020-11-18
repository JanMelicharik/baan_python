# ILUSTRACE GIBBSOVA VZORKOVACE
# Nastaveni cesty do domovske slozky
import sys
root_dir = "/".join(sys.argv[0].split("/")[:-2])
sys.path.insert(1, root_dir + "/support")

# importy pomocnych funkci
from progress_info import progress_bar

import numpy as np
from math import sqrt
from numpy.random import randn
import matplotlib.pyplot as plt

# ===== 1. Definovani parametru a pocatecni nastaveni =====
# Stredni hodnoty
mu = np.array([[0],[0]])
# Kovariancni matice
rho = 0.5
sigma = np.array([[1, rho], [rho, 1]])

s = 10_000      # pocet generovanych vzorku
s_1 = 5_000     # pocet ponechanych vzorku
s_0 = s - s_1   # pocet vyhozenych vzorku

# Matice pro generovane vzorky
theta = np.zeros((2, s+1))
theta[:,0] = np.array([-1000, 1000])

# ===== 2. Gibbsuv vzorkovac =====

for i in range(1,s+1):
    # generovani theta_1/theta_2 ~ N(mu_12,sigma_12)
    mu_12 = mu[0][0] + rho * (theta[1, i-1] - mu[1][0])
    sigma_12 = 1 - rho**2
    theta[0,i] = randn() * sqrt(sigma_12) + mu_12

    # generovani theta_2/theta_1 ~ N(mu_21,sigma_21)
    mu_21 = mu[1][0] + rho * (theta[0, i] - mu[0][0])
    sigma_21 = 1 - rho**2
    theta[1,i] = randn() * sqrt(sigma_21) + mu_21
    progress_bar(i, s+1)

# ===== 3. Grafy =====

fig, axs = plt.subplots(2, 2)

k = 10
axs[0, 0].plot(theta[0,0:k], theta[1,0:k], linewidth=0.5, marker="x")
axs[0, 0].set_title(f"Konvergence prvnich {k} kroku Gibbsova vzorkovace")
axs[0, 0].set_xlabel("theta_1")
axs[0, 0].set_ylabel("theta_2")

theta = theta[:,s_1:]
axs[0, 1].set_title(f"Sdruzena hustota na zaklade {s_1} vzorku")
axs[0, 1].scatter(theta[0], theta[1], marker=".", s=0.75)
axs[0, 1].set_xlabel("theta_1")
axs[0, 1].set_ylabel("theta_2")

axs[1, 0].set_title(f"Konvergence na zaklade 100. replikace")
axs[1, 0].plot(range(len(theta[0,::100])), theta[0,::100])
axs[1, 0].set_ylabel("theta_1")

axs[1, 1].set_title(f"Konvergence na zaklade 100. replikace")
axs[1, 1].plot(range(len(theta[1,::100])), theta[1,::100])
axs[1, 1].set_ylabel("theta_2")

# Vypis zakladnich statistik
means = np.mean(theta, axis=1)
covariances = np.cov(theta)
print("\nStredni hodnoty jednotlivych odhadu:")
for i in range(len(means)):
    print(f"- theta_{i+1}: {round(means[i], 4)}")
print(f"\nKovariancni matice:\n{np.matrix.round(covariances, 4)}")

fig.tight_layout()
plt.show()
