# Nastaveni cesty do domovske slozky
import sys
root_dir = "/".join(sys.argv[0].split("/")[:-2])
sys.path.insert(1, root_dir + "/support")

# importy pomocnych funkci
from norm_rnd import norm_rnd
from gamm_rnd_koop import gamm_rnd_koop
from progress_info import progress_bar
from geweke import geweke

# importy oficialnich knihoven
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from tabulate import tabulate
from scipy.io import loadmat
from numpy import transpose as t
from numpy.linalg import inv
from math import sqrt
from numpy import diag,\
                  ones,\
                  zeros


# Nacteni dat:
# time   ... cas cesty
# depart ... cas odchodu (v minutach od 6:30)
# reds   ... pocet cerevenych svetel na semaforu
# trains ... pocet vlaku, ktere musi nechat projek na Murrumbeena prejezdu

# Pro spravne nacteni datoveho souboru je treba spoustet script ze slozky 02_NLRM_s_nez...
# (popr. nastavit vlastni cestu)
raw_data = loadmat("data_commute.mat")

# Data jsou dictionary, kde hodnota je list hodnot
data = {}
for timeseries in ["reds", "time", "trains", "depart"]:
    data[timeseries] = raw_data[timeseries].transpose().tolist()[0]

# data jsou v matici v poradi: "reds", "time", "trains", "depart"
# data = np.array(data).transpose()

# ===== 1. Nastaveni apriornich hyperparametru a Gibbsova vzorkovace =====
# p(beta) ~ N(beta_0, v_0)
# p(h) ~ G(h_0, nu_0)
beta_0 = np.array([[30], [1], [3], [5]])
v_0 = diag([7.5**2, 0.25**2, 1**2, 2**2])
nu_0 = 40
s2_0 = 10**2
h_0 = 1 / s2_0

# Definice modelu
# Data musi byt zadefinovana jako matice i pokud se jedna o vektor
y = np.array([data["time"]]).transpose()
x = np.array([
        [1]*len(y),
        data["depart"],
        data["reds"],
        data["trains"]
    ]).transpose()

# Nastaveni Gibbsova vzorkovace
s = 50_000 + 1      # celkovy pocet generovanych vzorku + pocatecni hodnota
s_0 = 30_000 + 1    # pocet vyhozenych vzorku
s_1 = s - s_0       # pocet ponechanych vzorku

beta = np.array([[],[],[],[]])
beta = np.append(beta, beta_0, axis=1)

h = [h_0]

# ===== 2. Gibbsuv vzorkovac =====

nu_1 = len(y) + nu_0    # (4.9) dle Koop (2003)
for i in range(1,s):
    # 1. blok Gibbsova vzorkovace
    # podminena hustota p(beta|h,y) ~ N(beta_1, v_1)
    v_1 = inv(inv(v_0) + h[i-1] * (t(x) @ x))                   # (4.4) dle Koop (2003)
    beta_1 = v_1 @ (inv(v_0) @ beta_0 + h[i-1] * (t(x) @ y))    # (4.5)
    new_beta = beta_1 + norm_rnd(v_1)                           # (4.7)
    beta = np.append(beta, new_beta, axis=1)

    # 2. blok Gibbsova vzorkovace
    # podminena hustota p(h|beta,y) ~ G(h_1, nu_1)
    h_1 = inv((1 / nu_1 * (t(y - x @ beta[:,[i]]) @ (y - x @ beta[:,[i]]) + nu_0 * 1 / h_0)))[0][0]     # (4.10)
    new_h = gamm_rnd_koop(h_1, nu_1, 1)[0]    # (4.8)   (funkce gamm_rnd_koop vraci array - je treba rozbalit)
    h.append(new_h)
    progress_bar(i, s)

# ===== 3. Posteriorni analyza =====

# Vyhozeni prvnich s_0 vzorku
# Rozdilny zpusob mazani, protoze beta je matice a h je list
beta = np.delete(beta, range(s_0), axis=1)
h = h[s_0:]

# Graficke zobrazeni konvergence
k = 100     # delka kroku (neplest si s krokem ve funkci geweke)

fig, axs = plt.subplots(nrows=3, ncols=2)

for i in range(beta.shape[0]):
    row = 0 if i < 2 else 1
    col = i%2
    axs[row, col].plot(range(0, beta.shape[1], k), beta[i,::k])
    axs[row, col].set_title(f"beta_{i+1}")

axs[2,0].plot(range(0, len(h), k), h[::k])
axs[2,0].set_title("h")

axs[2,1].remove()

# Gewekova konvergencni diagnostika
cd_beta = geweke(beta)
cd_h = geweke(np.array([h]))

# ===== 4. Prezentace vysledku =====
# apriorni str. hodnoty a sm. odchylky
# beta_0, h_0 - apriorni stredni hodnoty
std_beta_0 = [sqrt(element) for element in np.diag(v_0)]
std_h_0 = sqrt(2 * h_0**2 / nu_0)

# posteriorni str. hodnoty a sm. odchylky
mean_beta_1 = np.mean(beta, axis=1)
mean_h_1 = np.mean(h)

std_beta_1 = np.std(beta, axis=1)
std_h_1 = np.std(h)

headers = ["Parametr", "Prior. mean", "Prior. std", "Post. mean", "Post. std", "CD"]
table = []

for i in range(beta_0.shape[0]):
    table.append(
        [
            f"beta_{i+1}",
            round(beta_0[i][0], 4),
            round(std_beta_0[i], 4),
            round(mean_beta_1[i], 4),
            round(std_beta_1[i], 4),
            round(cd_beta[i]["cd"], 4)
        ]
    )

table.append(
    [
        "h",
        round(h_0, 4),
        round(std_h_0, 4),
        round(mean_h_1, 4),
        round(std_h_1, 4),
        round(cd_h["cd"], 4)
    ]
)
print("\nVysledne statistiky:")
print(tabulate(table, headers, tablefmt="pretty"))

fig.tight_layout()
plt.show()