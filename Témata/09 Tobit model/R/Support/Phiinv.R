# Computes the standard normal quantile function of the vector x,
# 0<x<1.

Phiinv <- function(x) {
    val <- sqrt(2)*erfinv(2*x-1)
    return(val)
}
