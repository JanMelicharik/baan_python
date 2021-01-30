# Cviceni 04, priklad 01
rm(list = ls())
cat("\f")

# Definovani parametru a pocatecniho nastaveni G. vzorkovace
# stredni hodnoty
mu <- c(0,0)

# kovariancni matice
rho <- 0.5
Sigma <- matrix(c(1,rho,rho,1),nrow = 2,ncol = 2)

S <- 10000 # pocet generovanych vzorku
S1 <- 5000 # pocet ponechanych vzorku
S0 <- S - S1 # pocet vyhozenych vzorku

# matice pro generovane vzorky (matice samych nul)
theta <- matrix(0,ncol = 2,nrow = 2)
# pocatecni hodnota pro beh Gibsova vzorkovace = theta(0)
theta[,1] = c(-1000,1000)

# Gibbsuv vzorkovac
for(s in 2:(S+1)) {
    # generovani theta_1/theta_2 ~ N(mu_12,Sigma_12)
    mu_12 <- mu[1] + rho*(theta[2,s-1] - mu[2])
    Sigma_12 <- 1 - rho^2
    theta_12 <-  rnorm(1)*sqrt(Sigma_12) + mu_12
    # generovani theta_2/theta_1 ~ N(mu_21,Sigma_21)
    mu_21 <- mu[2] + rho*(theta_12 - mu[1])
    Sigma_21 <- 1 - rho^2
    theta_21 <- rnorm(1)*sqrt(Sigma_21) + mu_21
    theta <- cbind(theta[,1:(s-1)], c(theta_12,theta_21))
}

# zobrazime prvnich k kroku Gibbsova vzorkovace
k <- 10
par(mfrow = c(2,2))
plot(theta[1,1:k],theta[2,1:k],
     pch = 8,
     col = "blue",
     xlab = expression(~ theta[1]),
     ylab = expression(~ theta[2]))
lines(theta[1,1:k],theta[2,1:k],col = "blue")
title(c('Prvnich',as.character(k),'kroku Gibbsova vzorkovace'))

# Vyhodime prvnich S0+1 vzorku a vypocet statistik sdruzeneho rozdeleni
theta <- theta[,(S0+2):ncol(theta)]
theta_mean <- rowMeans(theta) # vypocet stredni hodnoty pres radky
theta_cov <- cov(t(theta))

cat('Vektor strednich hodnot\n')
print(theta_mean)

cat('\nKovariancni matice\n')
print(theta_cov)

# Graficke zobrazeni posteriorni hustoty
plot(theta[1,],theta[2,],
     pch = 20,
     col = "blue",
     xlab = expression(~ theta[1]),
     ylab = expression(~ theta[2]))
title(c('Sdruzena hustota na zaklade',as.character(S1),'vzorku'))

# graficke overeni konvergence (vykresleni kazdeho k-teho vzorku)
k <- 100
plot(theta[1,seq(1,ncol(theta),k)],
     type = "l",
     col = "blue",
     ylab = expression(~ theta[1]))
title(c('Konvergence na zaklade',as.character(k),'-te replikace'))

plot(theta[2,seq(1,ncol(theta),k)],
     type = "l",
     col = "blue",
     ylab = expression(~ theta[2]))
title(c('Konvergence na zaklade',as.character(k),'-te replikace'))



