# Gewekova konvergenční diagnostika
# vyuzivajici funkci momentg.m z
# LeSageho ekonometrickeho toolboxu
# theta ... rozmer S1 x k
# (S pocet vzorku, k pocet parametru)

Geweke <- function(theta) {
    if(is.vector(theta) == 1) {
        theta <- as.matrix(theta)
        if(nrow(theta) == 1) {
            theta <- t(theta)
        }
        S1 <- nrow(theta)
    } else if(is.vector(theta) != 1 & nrow(theta) == 1) {
        theta <- t(theta)
        S1 <- nrow(theta)
    } else {
        S1 <- nrow(theta)
    }

    smpl_A <- round(0.1*S1) # prvnich 10 % vzorku
    smpl_C <- round(0.6*S1) + 1 # poslednich 40 % vzorku

    # NSE pro vzorek A
    pom <- momentg(as.matrix(theta[1:smpl_A,]))
    mean_A <- pom$pmean
    nse_A <- pom$nse1

    # NSE pro vzorek B
    pom <- momentg(as.matrix(theta[(smpl_C+1):S1,]))
    mean_C <- pom$pmean
    nse_C <- pom$nse1

    # Gewekova konvergencni diagnostika
    CD <- (mean_A - mean_C)/(nse_A + nse_C)

    nse <- momentg(theta)

    results <- list()
    results$CD <- CD
    results$NSE <- nse$nse1

    # vysledek
    return(results)
}
