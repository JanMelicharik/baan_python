# PURPOSE: returns the inverse (quantile) at x of the chisq(n) distribution
# ---------------------------------------------------
#     USAGE: x = chis_inv(p,a)
#     where: p = a vector of probabilities
#            a = a scalar parameter
#     NOTE: chis_inv(x,n) = gamm_inv(p,a/2)*2
#    ---------------------------------------------------
#         RETURNS:
#            a vector x at each element of p from chisq(n) distribution
#
#                Anders Holtsberg, 18-11-93
#            Copyright (c) Anders Holtsberg

chis_inv <- function(p,a) {
    if(any(abs(2*p - 1) > 1)) {
        stop('chis_inv: a probability should be 0 <= p <= 1')
    }

    if(any(a <= 0)) {
        stop('chis_inv: dof is wrong')
    }
    res <- gamm_inv(p,a*0.5)*2
    return(res)
}
