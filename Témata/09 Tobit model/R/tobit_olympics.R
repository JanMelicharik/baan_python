# Tobit model - data o poctu ziskanych medaily (souhrnny panel)
rm(list = ls())
cat("\f")

# Nastaveni cest a nacteni dat
source("Support/norm_rnd.R")
source("Support/gamm_rnd_Koop.R")
source("Support/Geweke.R")
source("Support/momentg.R")
source("Support/stdn_inv.R")
source("Support/norm_inv.R")
source("Support/normt_rnd.R")
source("Support/normlt_rnd.R")
source("Support/normrt_rnd.R")
source("Support/stdn_cdf.R")
source("Support/norm_cdf.R")
library(pracma)

data <- read.csv("olympics.csv",sep = ';',dec = ',')

# Nacteni data
#  structure array data_olympics, 1610 pozorovani
#
#     country     country code
#     year        olympics year
#     gdp         gross domestic product, 1995 dollars
#     pop         population
#     gold        number of gold medals won
#     silver      number of silver medals won
#     bronze      number of bronze medals won
#     medaltot    total number of medals won
#     host        = 1 if host country
#     planned     = 1 if non-soviet planned
#     soviet      = 1 if soviet
#
# Data source: Andrew B. Bernard and Meghan R. Busse "Who wins the olympic games: Economic resources and medal totals,"
#     The Review of Economics and Statistics, February 2004, 86(1), 413-417
#
# Prevzato z Hill et al. (2007) doplneno o promennou:
#
#      share        = podil ziskanych medail v letech 1988, 1992 a 1996
#     (share = medaltot/738.*(year==88)+medaltot/815.*(year==92)+medaltot/842.*(year==96)

# Priprava dat
y_raw <- as.double(data$share)
X_raw <- cbind(rep(1,length(y_raw)),log(data$gdp),log(as.double(data$pop)))

# ocisteni o chyjejici hodnoty
ind_nona <- !is.na(y_raw)
y <- y_raw[ind_nona]
X <- X_raw[ind_nona,]
N <- length(y)
# Nastaveni apriornich hyperparametru a Gibbsova vzorkovace
# Apriornimi hyperparametry
# p(beta)~N(beta_0, V_0)
# p(h)~G(h_0,nu_0)
beta_0 <- c(0,0,0)
V_0 <- diag(c(0.2^2,0.01^2,0.01^2))
nu_0 <- 100
s2_0 <- 0.1^2
h_0 <- 1/s2_0

# Nastaveni Gibbsova vzorkovace
S <- 50000 + 1   # celkovy pocet generovanych vzorku + pocatecni hodnota
S_0 <- 30000 + 1 # pocet vyhozenych vzorku
S_1 <- S - S_0   # pocet ponechanych vzorku

beta <- mat.or.vec(length(beta_0),S)  # vzorky pro beta
h <- mat.or.vec(1,S)                  # vzorky pro h
y_ast <- mat.or.vec(N,S)      # vzorky pro y_ast (latentni data)

# nastaveni pocatecnich hodnot
beta[,1] <- beta_0
h[1] <- h_0
y_ast[,1] <- y


# Gibbsuv vzorkovac

# graficky ukazatel postupu simulace a zapnuti pocitadla casu
tic()

# Gibbsuv vzorkovac
cat('Gibbsuv vzorkovac - tobit model \n\n')
cat('0%      50%      100%\n')
cat('|')

# pomocna promenna pro graficke znazorneni prubehu
pom_step <- 0.05 # 0.05 = posun po 5%
pom_graph <- pom_step # pocitadlo prubehu simulace


for(s in 2:S) {
    # 1. blok Gibbsova vzorkovace
    # podminena hustota p(beta|h,y_ast) ~ N(beta_1,V_1)
    V_1 <- solve(solve(V_0) + h[s-1]*(t(X) %*% X)) # (4.4) dle Koop (2003)
    beta_1 <- V_1 %*% (solve(V_0) %*% beta_0 + h[s-1]*(t(X) %*% y_ast[,s-1]))
                                # (4.5) dle Koop (2003)
    beta[,s] <- beta_1 + norm_rnd(V_1) # (4.7) dle Koop (2003)

    # 2. blok Gibbsova vzorkovace
    # podminena hustota p(h|beta,y_ast) ~ G(h_1,nu_1)
    nu_1 <- N + nu_0    # (4.9)
    h_1 <- (1/nu_1*(t(y_ast[,s-1] - X %*% beta[,s]) %*% (
        y_ast[,s-1] - X %*% beta[,s]) + nu_0*1/h_0))^-1 # (4.10)

    h[s] <- gamm_rnd_Koop(h_1,nu_1,1)  # (4.8)

    # 3.blok Gibbsova vzorkovace
    # podminena hustota pro p(y_ast|beta,h,y)
    y_ast[,s] <- y
    # nalaezeni y == 0
    ind <- which(y == 0)
    # generovani nah. cisel z prava omezeneho normalniho rozdeleni
    c_r <- mat.or.vec(length(ind),1) # omezeni - nuly z prava
    y_ast[ind,s] <- normrt_rnd(X[ind,] %*% beta[,s],rep(1,length(ind))*1/h[s],c_r)

    # graficke zobrazeni prubehu po 5 %
    if(s/S >= pom_graph) {
        cat('|')
        pom_graph <- pom_graph + pom_step
    }
}
cat('\n')
cat('Done!\n')

# Posteriorni analyza
# vyhozeni prvnich S_0 vzorku
beta <- beta[,(S_0+1):S]
h <- h[(S_0+1):S]
y_ast <- y_ast[,(S_0+1):S]

# graficke zobrazeni konvergence
k <- 100 # delka kroku

par(mfrow = c(2,2))
for(ii in 1:length(beta_0)) {
    plot(beta[ii,seq(1,ncol(beta),k)],type = 'l',
         xlab = '',ylab = substitute(~ beta[ii], list(ii=ii)))
}
plot(h[seq(1,ncol(beta),k)],type = 'l',
     xlab = '',ylab = 'h')

# Gewekova konvergencni diagnostika
CD_beta <- Geweke(t(beta))
CD_h <- Geweke(t(h))


# Prezentace vysledku
# apriorni str. hodnodty a sm. odchylky
# beta_0, h_0 - apriorni stredni hodnoty
std_beta_0 <- sqrt(diag(V_0)) # vektor apriornich sm. odchylek
std_h_0 <- sqrt(2*h_0^2/nu_0) # apriorni sm. odchylka pro h

# posteriorni str. hodnoty a sm. dochylky
mean_beta_1 <- rowMeans(beta) # sloupcovy vektor radkovych prumeru
mean_h_1 <- mean(h) # vystup = skalarni velicina
std_beta_1 <- sqrt(rowMeans(beta^2) - mean_beta_1^2)
std_h_1 <- sqrt(mean(h^2) - mean_h_1)

# Vystup na obrazovku
cat('Parametr  prior m.   prior std.  post m.     post std.   CD\n')
cat('===========================================================\n')
for(ii in 1:length(beta_0)) {
    cat(sprintf('Beta %1.0f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f\n',
                ii,beta_0[ii],std_beta_0[ii],mean_beta_1[ii],
                std_beta_1[ii],CD_beta$CD[ii]))
}
cat(sprintf('h     %6.4f   %6.4f   %6.4f   %6.4f  6.4f\n',
            h_0,std_h_0,mean_h_1,std_h_1,CD_h$CD))


