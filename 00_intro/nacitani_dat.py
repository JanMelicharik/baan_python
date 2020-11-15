# Nacteni datoveho souboru capm2_data.mat.
# Priprava datoveho souboru capm2.mat.
# Vykresleni dat do grafu.

import pandas as pd
import numpy as np
import math

# Dve moznosti nacteni dat (u CSV souboru je nutne specifikovat
# oddelovat sloupcu a nastavit znak pro desetinnou carku).
excel_file = pd.read_excel("capm2_data.xlsx")
csv_file = pd.read_csv("capm2_data.csv", delimiter=";", decimal=",")
