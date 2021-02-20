# PURPOSE: computes the quantile (inverse of the CDF)
#           for each component of x with mean m, variance v
# ---------------------------------------------------
# USAGE: invp = norm_inv(x,m,v)
# where: x = variable vector (nx1)
#        m = mean vector (default=0)
#     v = variance vector (default=1)
# ---------------------------------------------------
# RETURNS: invp (nx1) vector
# ---------------------------------------------------
# SEE ALSO: norm_d, norm_rnd, norm_inv, norm_cdf
# ---------------------------------------------------
#
# Written by KH (Kurt.Hornik@ci.tuwien.ac.at) on Oct 26, 1994
# Copyright Dept of Probability Theory and Statistics TU Wien

norm_inv <- function(x,m = mat.or.vec(1,(nrow(x)*ncol(x))),v = t(rep(1,(nrow(x)*ncol(x))))) {
    r <- nrow(x)
    c <- ncol(x)

    s <- r*c

    x <- matrix(x,1,s)
    m <- matrix(m,1,s)
    v <- matrix(v,1,s)

    invp <- mat.or.vec(1,s)
    invp <- m + sqrt(v)*stdn_inv(x)
    invp <- matrix(invp,r,c)
    return(invp)
}
