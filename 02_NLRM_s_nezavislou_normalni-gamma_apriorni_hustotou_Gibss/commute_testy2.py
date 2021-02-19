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
from scipy.stats import norm
from scipy.stats import multivariate_normal as mvn
from numpy import transpose as t
from numpy.linalg import inv
from math import sqrt
from numpy import diag,\
                  ones,\
                  zeros

##### POZNAMKY:
# Predelat nacitani .mat souboru na CSV 


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

# Dalsi deklarace vektoru (pro Savage-Dickey pomer hustot)
sdd_nom = [0]                # jmenovatel pro SD pomer hustot (prvni prvek je nula, aby se predeslo chybe s indexovanim)
sde_nom = [0]                # citatel pro SD pomer hustot
count_post_restrict = 0      # pocitadlo nesplneni podminek apriornich restrikci na parametry

# ===== 2. Gibbsuv vzorkovac =====

print("Gibbsuv vzorkovac s apriorni restrikci:")
nu_1 = len(y) + nu_0    # (4.9) dle Koop (2003)
for i in range(1,s):
    # 1. blok Gibbsova vzorkovace
    # podminena hustota p(beta|h,y) ~ N(beta_1, v_1)
    v_1 = inv(inv(v_0) + h[i-1] * (t(x) @ x))                   # (4.4) dle Koop (2003)
    beta_1 = v_1 @ (inv(v_0) @ beta_0 + h[i-1] * (t(x) @ y))    # (4.5)

    check_restrict = True
    while check_restrict:
        new_beta = beta_1 + norm_rnd(v_1)                       # (4.7)
        count_post_restrict += 1
        if np.prod(new_beta >= 0):
            check_restrict = False
            count_post_restrict -= 1

    beta = np.append(beta, new_beta, axis=1)

    # 2. blok Gibbsova vzorkovace
    # Podminena hustota p(h|beta,y) ~ G(h_1, nu_1)
    h_1 = inv((1 / nu_1 * (t(y - x @ beta[:,[i]]) @ (y - x @ beta[:,[i]]) + nu_0 * 1 / h_0)))[0][0]     # (4.10)
    new_h = gamm_rnd_koop(h_1, nu_1, 1)[0]    # (4.8)   (funkce gamm_rnd_koop vraci array - je treba rozbalit)
    h.append(new_h)

    # Cast pro vypocet citatele SD pomeru hustot
    # d) beta_depart = 0 (druhy prvek vektoru beta)
    # p(b|h,y) ~ N(beta_1, v_1)
    v_1 = inv(inv(v_0) + h[i] * (t(x) @ x))                   # (4.4) dle Koop (2003)
    beta_1 = v_1 @ (inv(v_0) @ beta_0 + h[i] * (t(x) @ y))
    new_sdd_nom = norm.pdf(0, beta_1[1][0], sqrt(v_1[1][1]))
    sdd_nom.append(new_sdd_nom)

    # e) beta_reds = beta_trains
    # p(R*b|h,y) ~ N(R * beta_1, R * v_1 * R')
    r = np.array([[0, 0, 1, -1]])
    new_sde_nom = mvn.pdf(0, r @ beta_1, r @ v_1 @ t(r))
    sde_nom.append(new_sde_nom)

    progress_bar(i, s)

# ===== 3. Posteriorni analyza =====

# vyhozeni prvnich s_0 vzorku
beta = np.delete(beta, range(s_0), axis=1)
h = h[s_0:]
sdd_nom = sdd_nom[s_0:]
sde_nom = sde_nom[s_0:]

# Graficke zobrazeni konvergence
k = 100     # delka kroku (neplest si s krokem ve funkci geweke)
fig_1 = plt.figure()
axs = fig_1.subplots(nrows=3, ncols=2)

for i in range(beta.shape[0]):
    row = 0 if i < 2 else 1
    col = i%2
    axs[row, col].plot(range(0, beta.shape[1], k), beta[i,::k])
    axs[row, col].set_title(f"beta_{i+1}")

axs[2,0].plot(range(0, len(h), k), h[::k])
axs[2,0].set_title("h")

axs[2,1].remove()
fig_1.tight_layout()

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

# ===== 5. Test hypotezy ze beta_reds (beta_3) >= 2 =====
# pro ilustraci: zobrazeni rozdeleni beta_3 pomoci histogramu

fig_2 = plt.figure()
fig_2.suptitle("Posteriorni hustota parametru beta_3 (reds)")
axh = fig_2.add_subplot(1,1,1)
axh.hist(beta[2], bins=100, edgecolor="black", linewidth=1)
axh.set_xlabel("beta_3")

# pravdepodobnost beta_3 >= 2
pc = sum(beta[2] >= 2) / s_1
print("Pravdepodobnost beta_red >= 2 a odpovidajici Bayesuv faktor:")
print(f"- Pravdepodobnost = {round(pc, 4)}\n- Bayesuv faktor = {round(pc / (1 - pc), 4)}\n")

# ===== 6. Test hypotezy beta_depart (beta_2) == 0  =====
# Savage-Dickey pomer hustot
# Vyhodnotime p(beta_depart=0|M2)
# Apriorni hustota p(beta_depart)~N(beta_0(2),V_0(2,2))
sdd_denom = norm.pdf(0, beta_0[1][0], sqrt(v_0[1][1]))
s_sim = 50_000
count_prior_restrict = 0    # pocitadlo nesplneni apriornich omezeni

print("Simulace pro test hypotezy beta_depart (beta_2) == 0:")
for i in range(s_sim):
    check_restrict = True
    while check_restrict:
        pom = beta_0[1][0] + norm_rnd(np.array([[v_0[1][1]]]))
        count_prior_restrict += 1
        if pom >= 0:
            check_restrict = False
            count_prior_restrict -= 1
    progress_bar(i, s_sim)

# Integracni konstanta apriorni hustoty zohlednujici apriorni omezeni
int_prior_restrict = 1 / (1 - count_prior_restrict / (s_sim + count_prior_restrict))
# Korekce jmenovatele S-D pomeru hustot
sdd_denom = int_prior_restrict * sdd_denom
# Vypocet citatele na zaklade vystupu Gibbsova vzorkovace (rozsireneho)
e_sdd_nom = np.mean(sdd_nom)
# Integracni konstanta posteriorni hustoty zohlednujici apriorni omezeni
e_sdd_nom = 1 / (1 - count_post_restrict / (s + count_post_restrict - 1)) * e_sdd_nom
# Bayesuv faktor
sd_d = e_sdd_nom / sdd_denom

print("\nBayesuv faktor pro beta_depart = 0:")
print(f" - BF = {round(sd_d, 4)}")

# ===== 7. Test beta_reds (beta_3) == beta_trains (beta_4)  =====
# Linearni omezeni R*beta = r
# matice omezeni (na cely vektor parametru beta)
r = np.array([[0, 0, 1, -1]])
sde_denom = mvn.pdf(0, r @ beta_0, r @ v_0 @ t(r))
s_sim = 50_000
count_prior_restrict = 0

print("\nSimulace pro test beta_reds == beta_trains:")
for i in range(s_sim):
    check_restrict = True
    while check_restrict:
        pom = beta_0 + norm_rnd(v_0)
        count_prior_restrict += 1
        if np.prod(pom >= 0) == 1:
            check_restrict = False
            count_prior_restrict -= 1
    progress_bar(i, s_sim)

int_prior_restrict = 1 / (1 - count_prior_restrict / (s_sim + count_prior_restrict))
sde_denom = int_prior_restrict * sde_denom
e_sde_nom = np.mean(sde_nom)
e_sde_nom = 1 / (1 - count_post_restrict / (s + count_post_restrict - 1)) * e_sde_nom
sd_e = e_sde_nom / sde_denom

print("Bayesuv faktor pro beta_reds = beta_trains:")
print(f"- BF = {round(sd_e, 4)}")

# ===== 8. Vypocet intervalu nejvyssi posteriorni hodnoty =====

# vypocet kvantilu
hpdi_90_coeff = []
for coefficient in beta:
    hpdi_90_coeff.append(np.quantile(coefficient, [0.05,0.95]))

hpdi_90_h = np.quantile(h, [0.05,0.95])

headers = ["Parametr", "5% kvantil", "95% kvantil"]
table = []
for i in range(len(hpdi_90_coeff)):
    table.append([f"beta_{i+1}", round(hpdi_90_coeff[i][0], 4), round(hpdi_90_coeff[i][1], 4)])

table.append([f"h", round(hpdi_90_h[0], 4), round(hpdi_90_h[1], 4)])
print("\nTabulka 90% kvantilu pro odhadnute paramatry:")
print(tabulate(table, headers, tablefmt="pretty"))

# Vykresleni grafu
plt.show()
