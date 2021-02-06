import pdb
import numpy as np
from numpy import transpose as t
from numpy.linalg import inv

# Vytvareni jednoduchych matic: vstupem musi byt tuple nebo list nesouci rozmery dimenzi
# Vstupem do techto metod muze byt list, tuple, vektor ci matice (nezalezi na orientaci)
# Vytvoreni matice jednicek
matrix_ones = np.ones((4,3))        # 4 radky, 3 sloupce
# Vytvoreni matice nul
matrix_zeros = np.zeros([1,4])      # 1 radek, 4 sloupce - maticovy vektor (bezny vektor je v Pythonu neorientovany - 1 dimenze)
# Vytvoreni jednotkove matice
matrix_eye = np.diag([1]*4)         # Alternativa k funkci eye() v Matlabu
# Vytvoreni vlastni diagonalni matice
matrix_diag = np.diag((2,3,4,1))

# Vytvareni vlastnich matic: Doporucuje se pouzivat funkci np.array namisto np.matrix
# Matice musi byt zadana jako list listu, kazdy vnitrni list je radek
my_matrix = np.array(
    [
        [2,   3,   1,   4],
        [6,   2.1, 5,   0],
        [1.1, 4,   6,   5],
        [1,   4,   0,   2]
    ]
)
# Transpozice matice
my_t_matrix_1 = t(my_matrix)
my_t_matrix_2 = my_matrix.T
# Inverze matice: pouze ctvercove matice
my_inv_matrix = inv(my_matrix)

# Zobrazeni matice, poctu prvku a jejiho rozmeru
print(my_matrix)
print(f"Pocet prvku v my_matrix: {my_matrix.size}")
print(f"Dimenze my_matrix: {my_matrix.shape} (radky, sloupce)")

# Nasobeni matice skalarem
my_matrix_1 = my_matrix * 5
# Pricteni skalaru k matici
my_matrix_2 = my_matrix + 2

# Nasobeni matice vektorem
# Tato operace vraci vektor - NENI matice o rozmeru 1x4 nebo 4x1 - jedna se pouze o vektor (neorientovany)
my_vector_1 = my_matrix[1,:] @ my_matrix
# Nasobeni vektoru matici
my_vector_2 = my_matrix @ my_matrix[:,2]    

# Scitani matic
my_plus_matrix = my_matrix + my_t_matrix_1
# Odcitani matic
my_minus_matrix = my_matrix - my_t_matrix_1
# Nasobeni matic: NumPy zavedl pro nasobeni matic operator "@"
# Vysledkem je "temer" nulova matice - prvky mimo diagonalu jsou diky chybam
# behem vypoctu pocitace velmi nizke, ale ne zcela nulova (napr. 2.22e-16)
my_times_matrix = my_matrix @ my_inv_matrix

# Vyber jednotliveho prvku matice
my_scalar_1 = my_matrix[1][3]   # Vybira 2. radek a 4. sloupce (protoze Python indexuje od nuly)
my_scalar_2 = my_matrix[3][0]   # Vybira 4. radek a 1. sloupec

# Vyber radku matice -> Vektor
my_row_1 = my_matrix[2,:]
my_row_2 = my_matrix[0,:]
# Vyber sloupce matice -> Vektor
my_col_1 = my_matrix[:,3]
my_col_2 = my_matrix[:,0]

# Vyber radku matice -> Matice
row_matrix = my_matrix[[0],:]   # Vraci matici 1x4
# Vyber sloupce matice -> Matice
col_matrix = my_matrix[:,[2]]   # Vraci matici 4x1

# Nasobeni vektoru
# Skalarni soucin dvou vektoru: (v1, v2) * (u1, u2) = v1 * u1 + v2 * u2
my_scalar_3 = my_col_1 @ my_row_1
# Soucin vektoru: (v1, v2) * (u1, u2) = (v1 * u1, v2 * u2)
my_vector_3 = my_col_1 * my_row_1

# Nasobeni dvou vektoru tak, aby vznikla matice
# Zde neexistuje jednoduchy matematicky zapis a je ta treba pouzit metodu z knihovny NumPy - outer()
# Nelze tedy vzit dva vektory, jeden z nich transponovat a provest nasobeni tak, abychom dostali matici,
# protoze vektory jsou neorientovane.
my_new_matrix_1 = np.outer(my_row_1, my_col_1)
my_new_matrix_2 = np.outer(my_col_1, my_row_1)
# Kontrola, ze dve matice jsou stejne - pouhe porovnani matic vrati matici hodnot True/False.
# Tento vystup je treba agregovat pomoci funkce all().
same = (my_new_matrix_2.T == my_new_matrix_1).all()

# Vkladani prvku do matice
empty_matrix = np.zeros([4,4])
empty_matrix[0,0] = 5           # Prvek v prvnim radku a prvnim sloupci nahradi cislem 5
# Vkladani radku do matice
# Nasledujici radky produkuji stejny vysledek
empty_matrix[2,:] = [2,4,2,5]               # do radku matice vlozi list
empty_matrix[2,:] = (2,4,2,5)               # do radku matice vlozi tuple
empty_matrix[2,:] = np.array([2,4,2,5])     # do radku matice vlozi vektor
empty_matrix[2,:] = np.array([[2,4,2,5]])   # do radku matice vlozi maticovy vektor (matici 1x4) - nezalezi na orientaci

# Vkladani sloupce do matice
# Analogicky postup jako u vkladani radku
empty_matrix[:,2] = np.array([11,21,9,-12])

# Rozsirovani matic
# Tato operace ka matici muze priradit jakoukoliv matici, ktera splnuje shodu v rozmeru radku/sloupce.
# Pozor: Abychom matici rozsirili, je treba k ni pridavat radky/sloupce v maticovem zapisu - nelze do matice pridat vektor.
# Pozor: Zalezi na orientaci pridavaneho radku/sloupce - nelze pridat sloupec jako radek
# Pridani radku do matice
matrix_ones_1 = np.ones([4,4])
more_rows_1 = np.append(matrix_ones_1, [[2,2,2,2]], axis=0)     # axis=0 odpovida radkum
# Pozor: np.r_ neni funkce a tak nepouziva kulate zavorky "()", ale hranate "[]"
more_rows_2 = np.r_[matrix_ones_1, [[2,2,2,2]]]

# Pridani sloupce do matice
matrix_ones_2 = np.ones([4,4])
more_cols_1 = np.append(matrix_ones_2, [[3],[3],[3],[3]], axis=1)     # axis=1 odpovida sloupcum
more_cols_2 = np.c_[matrix_ones_2, [[3],[3],[3],[3]]]


# Pozastavi beh programu a umozni manipulaci a zobrazeni existujicich promennych
pdb.set_trace()     # Pro ukonceni programu staci do konzole "continue"
