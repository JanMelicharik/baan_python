# PURPOSE: returns the inverse of the cdf at p of the gamma(a) distribution
# --------------------------------------------------
#   USAGE: x = gamm_inv(p,a)
#   where: p = a vector of probabilities
#          a = a scalar parameter gamma(a)
#   --------------------------------------------------
#    RETURNS:
#         a vector x of the quantile at each element of p of the gamma(a) distribution
#   --------------------------------------------------
#
#     Anders Holtsberg, 18-11-93
#     Copyright (c) Anders Holtsberg

gamm_inv <- function(p,a) {
    if(any(abs(2*p - 1) > 1)) {
        stop('gamma_inv: a probability should be 0 <= p <= 1')
    }

    if(any(a <= 0)) {
        stop('gamma_inv: dof is wrong')
    }

    x <- max(a-1,0.1)
    dx <- 1
    while(any(abs(dx) > 256*(.Machine$double.eps)*max(x,1))) {
        dx <- (gamm_cdf(x,a) - p)/gamm_pdf(x,a)
        x <- x - dx
        x <- x + (dx - x)/2*(x < 0)
    }
    I0 <- which(p == 0)
    x[I0] <- 0
    I1 <- which(p == 1)
    x[I1] <- Inf
    return(x)
}
