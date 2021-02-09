# Gewekeho konvergencni diagnostika
# Autor: Bc. Jan Melicharik
# Zpracovano na zaklade Koop (2003), str. 66, rovnice (4.14)

# Tato funkce pouziva prepsanou funkci momentg.py

import sys
import numpy as np
from termcolor import colored
from math import sqrt

from support.momentg import momentg


def new_geweke(theta: np.ndarray, first: float = 0.1, last: float = 0.4, step: int = 10):
    """
        theta ... radkovy vektor ci matice (dvourozmerna) posteriornich vyberu
                  Gibbsova vzorkovace
    """
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
    
    s_1 = dimensions[0]
    smpl_a = round(0.1 * s_1)
    smpl_c = round(0.6 * s_1)

    pom = momentg(theta[:smpl_a, :])
    mean_a = pom["pmean"]
    nse_a = pom["nse1"]

    pom = momentg(theta[smpl_c:, :])
    mean_c = pom["pmean"]
    nse_c = pom["nse1"]

    cd = (mean_a - mean_c) / (nse_a + nse_c)

    pom = momentg(theta)
    nse = pom["nse1"]

    return {
        "cd": cd,
        "nse": nse
    }
