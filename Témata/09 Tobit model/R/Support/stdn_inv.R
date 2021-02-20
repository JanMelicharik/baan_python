# PURPOSE: computes the quantile (inverse of the CDF)
#          for each component of x with mean 0, variance 1
# ---------------------------------------------------
# USAGE: ninv = stdn_inv(x)
# where: x = variable vector (nx1)
# ---------------------------------------------------
# RETURNS: ninv = (nx1) vector containing quantiles at each x-element
# ---------------------------------------------------
#
# Written by KH (Kurt.Hornik@ci.tuwien.ac.at)

stdn_inv <- function(x) {
    ninv <- sqrt(2)*erfinv(2*x-1)
    return(ninv)
}




