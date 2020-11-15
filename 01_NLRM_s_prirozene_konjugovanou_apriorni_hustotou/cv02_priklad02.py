# Nastaveni cesty do domovske slozky
import sys
root_dir = "/".join(sys.argv[0].split("/")[:-2])
sys.path.insert(1, root_dir + "/support")

from numpy.random import normal,\
                         uniform,\
                         gamma
from math import sqrt
from statistics import mean
from tabulate import tabulate

from gamm_rnd_koop import gamm_rnd_koop

import matplotlib.pyplot as plt

headers = ["Velikost vzorku", "E", "D","std.err."] # zahlavi tabulek

# Definovani vektoru s pocty generovanych vzorku
MC = [[10, 3], [100, 10], [10_000, 50]]

# a) MC simulace z N(1,4)
fig_1, axs_1 = plt.subplots(3) # nastaveni pro vykreslovani
i = 0

# Nastaveni tabulky
table_1 = [headers]

for sample_size, bins in MC:
    sample = normal(loc=1, scale=sqrt(4), size=sample_size)
    E = mean(sample)                # vyberovy prumer jako odhad str. hodnoty
    D = mean(sample ** 2) - E ** 2  # vyberovy rozptyl

    table_1.append([sample_size, E, D, sqrt(D)])    # doplneni radku do tabulky

    axs_1[i].hist(
        sample,
        bins,
        facecolor="blue", 
        alpha=0.5,
        edgecolor="black",
        linewidth=1)
    axs_1[i].set_title(f"MC = {sample_size}")
    fig_1.suptitle("Vzorky z N(1,4)", fontsize=16)
    
    i += 1

print(f"\nSimulace rozdeleni N(1,4):\n{'='*35}")
print(tabulate(table_1, headers="firstrow"))

# b) MC simulace z U(2,5)
fig_2, axs_2 = plt.subplots(len(MC))
i = 0

table_2 = [headers]

for sample_size, bins in MC:
    sample = uniform(low=2, high=5, size=sample_size)
    E = mean(sample)
    D = mean(sample ** 2) - E ** 2

    table_2.append([sample_size, E, D, sqrt(D)])

    axs_2[i].hist(
        sample,
        bins,
        facecolor="blue", 
        alpha=0.5,
        edgecolor="black",
        linewidth=1)
    axs_2[i].set_title(f"MC = {sample_size}")
    fig_2.suptitle("Vzorky z U(2,5)", fontsize=16)

    i += 1

print(f"\nSimulace rozdeleni U(2,5):\n{'='*35}")
print(tabulate(table_2, headers="firstrow"))

# c) MC simulace z G(2,10) - podle Koopa
fig_3, axs_3 = plt.subplots(len(MC))
i = 0

table_3 = [headers]

for sample_size, bins in MC:
    sample = gamm_rnd_koop(2, 10, sample_size)
    E = mean(sample)
    D = mean(sample ** 2) - E ** 2

    table_3.append([sample_size, E, D, sqrt(D)])

    axs_3[i].hist(
        sample,
        bins,
        facecolor="blue", 
        alpha=0.5,
        edgecolor="black",
        linewidth=1)
    axs_3[i].set_title(f"MC = {sample_size}")
    fig_3.suptitle("Vzorky z G(2,10)", fontsize=16)

    i += 1

print(f"\nSimulace rozdeleni G(2,10):\n{'='*35}")
print(tabulate(table_3, headers="firstrow"))

# Zaverecne vykresleni

fig_1.tight_layout()
fig_2.tight_layout()
fig_3.tight_layout()
plt.show()
