# Nastaveni cesty do domovske slozky
import sys
root_dir = "/".join(sys.argv[0].split("/")[:-2])
sys.path.insert(1, root_dir + "/support")

# Importy ze slozky support
from progress_info import progress_bar

# Obecne importy
import time

from numpy.random import uniform
from tabulate import tabulate

start = time.time()

# ===== Hra Alice a Boba =====

## 0 - Nastaveni parametru pro ruzne specifikace hry
n = 8       # pocet her
a_win = 5   # pocet vyher Alice (musi platit: a_win <= n)

## 1 - Deklarace promennych
S = 1_000_000   # pocet simulaci
r = 0           # pocitadlo stavu n/a_win
b_win = 0       # stredni hodnota vyhry Boba pri stavu n/a_win

step = round(S/100)

## 2 - Simulace hry
for i in range(S):
    p = uniform()               # pocatecni rozdeleni stolu p~U(0,1)
    y = uniform(size=n)         # vektor rozmeru 8x1 s nezavislymi tahy z U(0,1)
    if sum(y < p) == a_win:
        r += 1                  # pocitadlo stavu n/a_win inkrementujeme o 1
        b_win += (1 - p)**3     # scitani pravdepodobnosti vyhry Boba
                                # (z predchozich behu pri stavu n/a_win)
    progress_bar(i, S)

## 3 - Vypocet a prezentace vysledku
b_win_expected = b_win / r          # Stredni (ocekavana) hodnota pravd. vyhry Boba
                                    # pri stavu n/a_win
a_win_expected = 1 - b_win_expected # Stredni (ocekavana) hodnota pravd. vyhry Alice
                                    # pri stavu n/a_win

print(f"\nPrumerna pravd. vyhry jednotlivych hracu pri stavu {a_win}:{n-a_win}:")
print(tabulate(
    [
        ["Jmeno", "Pravd."],
        ["Alice", a_win_expected],
        ["Bob",b_win_expected]
    ], headers="firstrow"))

print(f"\nPocet stavu {a_win}:{n-a_win}: {r}")
print(f"\nRelativni zastoupeni poctu stavu {a_win}:{n-a_win} na celkovem poctu simulaci: {r/S}")

end=time.time()
print(f"\nProces trval {round(end - start, 2)}s.")
