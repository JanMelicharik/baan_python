# Generator nahodnych cisel z G(mu,nu) rozdeleni
# dle Koop (2003) - vyuziva gamm_rnd z Econ. toolboxu
# mu ... stredni hodnota
# nu ... pocet stupnu volnosti
# m x 1 ... rozmer nahodneho vektoru y

gamm_rnd_Koop <- function(mu,nu,m){
    A <-  nu/2
    B <-  2*mu/nu
    res <- B*rgamma(m,A)

    return(res)
}
