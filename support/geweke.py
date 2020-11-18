# Gewekeho konvergencni diagnostika
# Autor: Bc. Jan Melicharik
# Zpracovano na zaklade Koop (2003), str. 66, rovnice (4.14)

# Zavislost vzorku je "vyresena" omezenym vyberem z jednotlivych casti vzorku.
# Vstupni vzorek se rozdeli na dve casti na koncich vzorku podle vstupu first a last.
# Z kazde teto casti se pro vypocet bere pouze zlomek odhadu (podle vstupu "step").

import sys
import numpy as np
from math import sqrt


def geweke(theta: np.ndarray, first: float = 0.1, last: float = 0.4, step: int = 10):
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

    if max(dimensions) * min(first, last) / step < 20:
        # Omezeni aby slo provest odhady na vektorech o delce alespon 20
        print("Nedostatecny pocet pozorovani pro dane nastaveni kroku.")

    number_of_coefficients = np.min(dimensions)
    is_column_wise = dimensions.index(number_of_coefficients)
    if is_column_wise:
        # Pokud je vstupem matice, kde koeficienty jsou jednotlive sloupce, transponuje ji
        theta = np.transpose(theta)

    if 1 not in dimensions:
        # provede funkci geweke po jednotlivych koeficientech
        number_of_coefficients = np.min(dimensions)
        return [geweke(theta[[coefficient]], first, last) for coefficient in range(number_of_coefficients)]

    n_obs = theta.shape[1]

    s_a = round(n_obs * first)               # prvnich 10%
    s_c = n_obs - round(n_obs * last)        # poslednich 40%

    a = theta[0][:s_a:step]
    c = theta[0][s_c::step]

    e_a = np.mean(a)
    e_c = np.mean(c)

    n_a = len(a)
    n_c = len(c)

    nse_a = np.std(a) / sqrt(n_a)
    nse_c = np.std(c) / sqrt(n_c)

    cd = (e_a - e_c) / (nse_a + nse_c)

    results = {
        "theta": theta,
        "s_a":   s_a,
        "s_c":   s_c,
        "a":     a,
        "c":     c,
        "n_a":   n_a,
        "n_c":   n_c,
        "e_a":   e_a,
        "e_c":   e_c,
        "nse_a": nse_a,
        "nse_c": nse_c,
        "cd":    cd
    }

    return results