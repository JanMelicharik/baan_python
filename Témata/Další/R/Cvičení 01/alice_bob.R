# Muj prvni skript
rm(list = ls())
cat("\f") # prikaz vycisti command window

# 1. Hra Alice a Boba - deklarace promennych a parametru
S <- 1000000 # pocet opakovani cele hry
r <- 0 # pocitadlo stavu 5:3
E_Bwin <- 0 # stredni hodnota vyhry Boba pri stavu 5:3

# 2. Simulace hry
for (s in 1:S) {
    p <- runif(1) # pocatecni rozdeleni stolu
         # generujeme p ~ U(0,1)
    # zahrani 8 her
    y <- runif(8) # vektor rozmeru 8x1 s nezavislymi U(0,1)
    # overime stav 5:3 => Alice musi vyhrat prave 5x
    # y < p ... vrati vektor 0/1, kde 1 = vyhrala Alice
    # sum(y<p) ... spocita pocet vyher Alice -> porovname
    #              s 5 sum ... secte prvky vektoru
    if (sum(y < p) == 5) {
        r <- r + 1 # pocitadlo stavu 5:3 se posune o jedna
        # scitame pravdepodobnosti vyhry Boba
        # (z predchoziho behu pri stavu 5:3)
        E_Bwin <- E_Bwin + (1 - p)^3
    }
}


# Vypocet a prezentace pozadovanych statistik
# Stredni (ocekavana) hodnota pravdepodobnosti vyhry Boba
# pri stavu 5:3
E_Bwin <- E_Bwin/r
# Stredni (ocekavana) hodnota pravdepodobnosti vyhry Alice
# pri stavu 5:3
E_Awin <- 1 - E_Bwin

# Jednoduchy vypis vysledku na obrazovku
print('Prumerna pravdepodobnost vyhry Boba pri stavu 5:3')
print(E_Bwin)
print('Prumerna pravdepodobnost vyhry Alice pri stavu 5:3')
print(E_Awin)
print('Ferovy podil sanci A:B pro rozdeleni vyhry za stavu 5:3')
print(E_Awin/E_Bwin)
print('Pocet stavu 5:3')
print(r)
print('Relativni zastoupeni poctu stavu 5:3 na celkovem poctu simulaci')
print(r/S)

# Graficke zobrazeni vysledku
plot(1, E_Bwin, xlim = range(c(0,3)), ylim = range(c(0,1)), pch = 8, col = "red")
points(2, E_Awin, pch = 3, col = "green")

