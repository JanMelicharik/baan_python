# Muj prvni skript
rm(list = ls())
cat("\f") # prikaz vycisti command window


# 1. Hra Alice a Boba - deklarace promennych a parametru
S <- 1000000 # pocet opakovani cele hry
r <- 0 # pocitadlo stavu a:b
E_Bwin <- 0 # stredni hodnota vyhry Boba pri stavu a:b
k <- 29 # celkovy pocet her
a <- 19 # pocet vyher Alice
b <- 10 # pocet vyher Boba


# 2. Simulace hry
for (s in 1:S) {
    p <- runif(1) # pocatecni rozdeleni stolu
    # generujeme p ~ U(0,1)
    # zahrani k her
    y <- runif(k) # vektor rozmeru kx1 s nezavislymi U(0,1)
    # overime stav a:b => Alice musi vyhrat prave ax
    # y < p ... vrati vektor 0/1, kde 1 = vyhrala Alice
    # sum(y<p) ... spocita pocet vyher Alice -> porovname
    #              s a sum ... secte prvky vektoru
    if (sum(y < p) == a) {
        r <- r + 1 # pocitadlo stavu a:b se posune o jedna
        # scitame pravdepodobnosti vyhry Boba
        # (z predchoziho behu pri stavu a:b)
        E_Bwin <- E_Bwin + (1 - p)^b
    }
}


# Vypocet a prezentace pozadovanych statistik
# Stredni (ocekavana) hodnota pravdepodobnosti vyhry Boba
# pri stavu a:b
E_Bwin <- E_Bwin/r
# Stredni (ocekavana) hodnota pravdepodobnosti vyhry Alice
# pri stavu a:b
E_Awin <- 1 - E_Bwin

# Jednoduchy vypis vysledku na obrazovku
cat('Prumerna pravdepodobnost vyhry Boba pri stavu ',a,':',b,'\n')
print(E_Bwin)
cat('Prumerna pravdepodobnost vyhry Alice pri stavu ',a,':',b,'\n')
print(E_Awin)
cat('Ferovy podil sanci A:B pro rozdeleni vyhry za stavu ',a,':',b,'\n')
print(E_Awin/E_Bwin)
cat('Pocet stavu ',a,':',b,'\n')
print(r)
cat('Relativni zastoupeni poctu stavu ',a,':',b,' na celkovem poctu simulaci\n')
print(r/S)

# Graficke zobrazeni vysledku
plot(1, E_Bwin, xlim = range(c(0,3)), ylim = range(c(0,1)), pch = 8, col = "red")
points(2, E_Awin, pch = 3, col = "green")

