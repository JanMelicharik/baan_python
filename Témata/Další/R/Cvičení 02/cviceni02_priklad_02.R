# Cviceni 02, priklad 02
rm(list = ls())
cat("\f")

# Definovani vektoru s pocty generovanych vzorku
MC <- c(10,100,10000)

# a) MC simulace z N(1,4)
cat('Simulace rozdeleni N(1,4)\n')
cat('=========================\n')

par(mfrow = c(3,1)) # nastavi plochu pro vykreslovani grafu
                    # pro vykresleni vice grafu
                    # aby byly vykreslovany grafy po jednom
                    # je nutne nastavit mfrow = c(1,1)
for(s in 1:length(MC)) {
    # generovani vzorku s N(1,4)
    smpl <- 1 + 2*rnorm(MC[s]) # rozptyl = 4
                               # sm. odchylka = 2
    # pocitani statistik
    E <- mean(smpl) # vyberovy prumer jako odhad
                    # stredni hodnoty
    D <- mean(smpl^2) - E^2
    # vystup na obrazovku (formatovany)
    cat('MC =',MC[s],'E =',E,'D =',D,'sd =',sqrt(D),'\n')

    # vykresleni pomoci histogramu
    hist(smpl)
    title(c('MC = ',as.character(MC[s])))
}

# b) MC simulace z U(2,5)
cat('\n\n')
cat('Simulace rozdeleni U(2,5)\n')
cat('=========================\n')

for(s in 1:length(MC)) {
    # generovani vzorku s U(2,5)
    smpl <- 2 + (5-2)*runif(MC[s])
    # pocitani statistik
    E <- mean(smpl) # vyberovy prumer jako odhad
    # stredni hodnoty
    D <- mean(smpl^2) - E^2
    # vystup na obrazovku (formatovany)
    cat('MC =',MC[s],'E =',E,'D =',D,'sd =',sqrt(D),'\n')

    # vykresleni pomoci histogramu
    hist(smpl)
    title(c('MC = ',as.character(MC[s])))
}

# c) MC simulace z G(2,10) ... dle Koopa
cat('\n\n')
cat('Simulace rozdeleni G(2,10)\n')
cat('=========================\n')

source("gamrnd_Koop.R")
for(s in 1:length(MC)) {
    # generovani vzorku s U(2,5)
    smpl <- gamrnd_Koop(2,10,MC[s])
    # pocitani statistik
    E <- mean(smpl) # vyberovy prumer jako odhad
    # stredni hodnoty
    D <- mean(smpl^2) - E^2
    # vystup na obrazovku (formatovany)
    cat('MC =',MC[s],'E =',E,'D =',D,'sd =',sqrt(D),'\n')

    # vykresleni pomoci histogramu
    hist(smpl)
    title(c('MC = ',as.character(MC[s])))
}
