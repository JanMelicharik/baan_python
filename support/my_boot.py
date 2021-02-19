# Vazeny bootstrap
#  index ... S*1 vektor bootstrapovych indexu
#  y ....... N x 1 vektor bootstrapovaneho vzorku
#  S ....... velikost bootstrapoveho vyberu
#  w ....... bootstrapove vahy (optional)

from support.progress_info import progress_bar
import numpy as np

def my_boot(y: np.ndarray, s: int, w=None):
    n = y.size
    # if not w:
    #     w = 1/n
    # Vektor indexu
    boot_index = [0] * s
    # Vytvoreni bodu intervalu, kde sirka dilcich intervalu reprezentuje stanovenou vahu
    w_int = np.cumsum(w)
    # Nastaveni seedu (muze byt libovolne cislo) - zakomentovano (neni nezbytne nutne)
    # np.random.seed(420)
    print("\nBootstrapping:")
    for i in range(s):
        threshold = np.random.rand() * w_int[-1]
        boot_index[i] = int(np.sum(w_int < threshold))   # konverze z float na int protoze se pozdeji pouziva jako index
        progress_bar(i, s)

    return boot_index
