# Matice v Pythonu
NumPy používá tzv. **broadcasting** kdyz pracuje s maticemi (`numpy.array()`), které nejsou stejného rozměru. V případě adresování jednotlivých částí matice (sloupce, řádky) je tak jištější adresovat sloupec/řádek jeho indexem v hranatých závorkách; např. (`my_matrix[:,[2]]`). [Více informací o broadcastingu](https://numpy.org/devdocs/user/theory.broadcasting.html).
