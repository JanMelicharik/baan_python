# Vraci vektor M x 1 nezavislych vyberu
# z gama rozdeleni dle znaceni Koop(2003)
# mu ... stredni hodnota
# nu ... pocet stupnu volnosti

gamrnd_Koop <- function(mu,nu,m){
    A <-  nu/2
    B <-  2*mu/nu
    res <- B*rgamma(m, shape = A, scale = B)

    return(res)
}
