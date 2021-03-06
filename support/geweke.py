# Gewekeho konvergencni diagnostika
# Autor: Bc. Jan Melicharik
# Zpracovano na zaklade Koop (2003), str. 66, rovnice (4.14)

# Zavislost vzorku je "vyresena" omezenym vyberem z jednotlivych casti vzorku.
# Vstupni vzorek se rozdeli na dve casti na koncich vzorku podle vstupu first a last.
# Z kazde teto casti se pro vypocet bere pouze zlomek odhadu (podle vstupu "step").

import sys
import numpy as np
from termcolor import colored
from math import sqrt


def geweke(theta: np.ndarray, first: float = 0.1, last: float = 0.4, step: int = 10):
    if first <= 0 or first >= 1:
        print("Pomer prvni casti vzorku musi byt z intervalu (0,1).")
        sys.exit(1)

    if last <= 0 or last >= 1:
        print("Pomer posledni casti vzorku musi byt z intervalu (0,1).")
        sys.exit(1)

    if first + last > 1:
        print("Soucet pomeru prvni a posledni casti vzorku nesmi byt vetsi nez 1. Napr. prvni: 0.1, posledni: 0.4")
        sys.exit(1)

    if step < 5 and step > 0:
        print("Pri nizkych hodnotach kroku mohou byt vysledky zkresleny vzajemnou korelaci.")

    if not isinstance(step, int) or step < 1:
        print("Krok musi byt cele kladne cislo. Napr. step = 10")
        sys.exit(1)

    dimensions = theta.shape
    if len(dimensions) != 2:
        print("Vlozeny objekt neni dvourozmerna matice. I vektory musi byt takto zadany.")
        sys.exit(1)

    if dimensions[0] < dimensions[1]:
        print(colored("Vstupem musi byt matice, jejiz sloupce jsou jednotlive promenne.", "red"))
        sys.exit(1)

    n_obs = theta.shape[0]
    a = round(first * n_obs)
    c = round((1 - last) * n_obs)

    smpl_a = theta[:a:step, :]
    smpl_c = theta[c::step, :]

    e_a = np.mean(smpl_a, axis=0)
    e_c = np.mean(smpl_c, axis=0)

    n_a = smpl_a.shape[0]
    n_c = smpl_c.shape[0]

    nse_a = np.std(smpl_a, axis=0) / sqrt(n_a)
    nse_c = np.std(smpl_c, axis=0) / sqrt(n_c)

    nse = np.std(theta, axis=0) / sqrt(n_obs)

    cd = (e_a - e_c) / (nse_a + nse_c)

    results = {
        "theta":  theta,
        "smpl_a": smpl_a,
        "smpl_c": smpl_c,
        "n_a":    n_a,
        "n_c":    n_c,
        "e_a":    e_a,
        "e_c":    e_c,
        "nse_a":  nse_a,
        "nse_c":  nse_c,
        "nse":    nse,
        "cd":     cd
    }

    return results
