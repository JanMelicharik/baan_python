# Funkce pro bayesianskou analyzu NLRM (prirozene konjugovana apriorni hustota)
# y ... vektory vysvetlujici promenne N x 1
# X ... matice plánu N x k
# beta_0, V_0, h_0, nu_0 ... apriorní hyperparametry
#                            z p(beta,h)~NG(beta_0,V_0,h_0,nu_0)
# beta_1, V_1, h_1, nu_1 ... aposteriorní hyperparametry
#                            z p(beta,h|y)~NG(beta_0,V_0,h_0,nu_0)
# b0_cov, b0_std, h0_std ... apriorní kovarianční matice
#			     a vektory apriorních směrodatných odchylek parametrů
# b1_cov, b1_std, h1_std ... posteriorní kovarianční matice
#			     a vektory posteriorních směrodatných odchylek parametrů
# log_ML 		 ... logaritmus marginální věrohodnosti modelu

my_NLRM <- function(y,X,beta_0,V_0,h_0,nu_0){
    if(is.matrix(X) == TRUE){
        # Pocitani charakteristik apriornich hustot
        b0_cov <- nu_0*h_0^(-1)/(nu_0 - 2)*V_0 # analogie (3.16) z Koop (2003)
        b0_std <- sqrt(diag(b0_cov))
        # apriorni sm. odchylka pro h
        h0_std <- sqrt(2*h_0/nu_0)

        # Pocitani aposteriornich hyperparametru
        N <- length(y) # pocet pozorovani
        # odhady OLS
        b_OLS <- (solve(t(X) %*% X) %*% t(X) %*% y) # (3.5) z Koop (2003)
        nu_OLS <- N - ncol(X) # (3.4) z Koop (2003)
        s2_OLS <- as.vector(1/nu_OLS*t(y - X %*% b_OLS) %*% (y - X %*% b_OLS)) # (3.6) z Koop (2003)

        # Aposteriorni hyperparametry
        V_1 <- solve(solve(V_0) + t(X) %*% X) # (3.10) z Koop (2003)
        beta_1 <- V_1 %*% (solve(V_0) %*% beta_0 + t(X) %*% X %*% b_OLS) # (3.11) z Koop (2003)
        nu_1 <- nu_0 + N # (3.12) z Koop (2003)
        h_1 <- as.vector(nu_1*(nu_0*1/h_0 + as.vector(nu_OLS*s2_OLS) + t(b_OLS - beta_0) %*% solve(V_0 + solve(t(X) %*% X)) %*% (b_OLS - beta_0))^(-1))
        # (3.13) z Koop (2003), kdy h_1 = s2_1^-1

        # Aposteriorni kovariancni matice a smerodatne odchylky
        b1_cov <- nu_1*h_1^(-1)/(nu_1 - 2)*V_1 # (3.16) z Koop (2003)
        b1_std <- sqrt(diag(b1_cov))

        # Aposteriorni smerodatna odchylka pro presnost chyby h
        h1_std <- sqrt(2*h_1/nu_1) # odmocnina z (3.19) z Koop (2003)

        # Logaritmus marginalni verohodnosti
        # log konstanty a marginalni verohodnosti
        log_c <- lgamma(nu_1/2) + nu_0/2*log(nu_0/h_0) - lgamma(nu_0/2) - N/2*log(pi) # logaritmus (3.35) z Koop (2003)
        log_ML <- log_c + 1/2*(log(det(V_1)) - log(det(V_0))) - nu_1/2*log(nu_1/h_1) # logaritmus (3.34) z Koop (2003)
    }
    else {
        # Pocitani charakteristik apriornich hustot
        b0_cov <- nu_0*h_0^(-1)/(nu_0 - 2)*V_0 # analogie (3.16) z Koop (2003)
        b0_std <- sqrt(b0_cov)
        # apriorni sm. odchylka pro h
        h0_std <- sqrt(2*h_0/nu_0)

        # Pocitani aposteriornich hyperparametru
        N <- length(y) # pocet pozorovani
        # odhady OLS
        b_OLS <- as.vector(solve(t(X) %*% X) %*% t(X) %*% y) # (3.5) z Koop (2003)
        nu_OLS <- N - 1 # (3.4) z Koop (2003)
        s2_OLS <- as.vector(1/nu_OLS*t(y - X*b_OLS) %*% (y - X*b_OLS)) # (3.6) z Koop (2003)

        # Aposteriorni hyperparametry
        V_1 <- as.vector(solve(solve(V_0) + t(X) %*% X)) # (3.10) z Koop (2003)
        beta_1 <- V_1*as.vector(V_0^(-1)*beta_0 + as.vector((t(X) %*% X)*b_OLS)) # (3.11) z Koop (2003)
        nu_1 <- nu_0 + N # (3.12) z Koop (2003)
        h_1 <- as.vector(nu_1*(nu_0*1/h_0 + as.vector(nu_OLS*s2_OLS) + t(b_OLS - beta_0) %*% solve(V_0 + solve(t(X) %*% X))*(b_OLS - beta_0))^(-1))
        # (3.13) z Koop (2003), kdy h_1 = s2_1^-1

        # Aposteriorni kovariancni matice a smerodatne odchylky
        b1_cov <- nu_1*h_1^(-1)/(nu_1 - 2)*V_1 # (3.16) z Koop (2003)
        b1_std <- sqrt(b1_cov)

        # Aposteriorni smerodatna odchylka pro presnost chyby h
        h1_std <- sqrt(2*h_1/nu_1) # odmocnina z (3.19) z Koop (2003)

        # Logaritmus marginalni verohodnosti
        # log konstanty a marginalni verohodnosti
        log_c <- lgamma(nu_1/2) + nu_0/2*log(nu_0/h_0) - lgamma(nu_0/2) - N/2*log(pi) # logaritmus (3.35) z Koop (2003)
        log_ML <- log_c + 1/2*(log(V_1) - log(V_0)) - nu_1/2*log(nu_1/h_1) # logaritmus (3.34) z Koop (2003)
    }
    # Vystup funkce do listu res
    res <- list(
        # puvodni data
        y = y,
        X = X,

        # apriorni hyperparametry
        beta_0 = beta_0,
        h_0 = h_0,
        V_0 = V_0,
        nu_0 = nu_0,

        b0_cov = b0_cov,
        b0_std = b0_std,
        h0_std = h0_std,

        # aposteriorni hyperparametry a kovariance
        beta_1 = beta_1,
        h_1 = h_1,
        V_1 = V_1,
        nu_1 = nu_1,

        b1_cov = b1_cov,
        b1_std = b1_std,
        h1_std = h1_std,

        # logaritmus maximalni verohodnosti
        log_ML = log_ML
    )

    return(res)
}
