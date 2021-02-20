# PURPOSE: computes the cumulative normal distribution
#           for each component of x with mean m, variance v
# ---------------------------------------------------
# USAGE: cdf = norm_cdf(x,m,v)
# where: x = variable vector (nx1)
#        m = mean vector (default=0)
#        v = variance vector (default=1)
# ---------------------------------------------------
# RETURNS: cdf (nx1) vector
# ---------------------------------------------------
#
# Written by TT (Teresa.Twaroch@ci.tuwien.ac.at) on Jun 3, 1993
# Updated by KH (Kurt.Hornik@ci.tuwien.ac.at) on Oct 26, 1994
# Copyright Dept of Probability Theory and Statistics TU Wien
# Updated by James P. Lesage, jpl@jpl.econ.utoledo.edu 1/7/97


norm_cdf <- function(x,m = mat.or.vec(nrow(as.matrix(x)),1),v = rep(1,nrow(as.matrix(x)))) {
    r <- nrow(as.matrix(x))
    c <- ncol(as.matrix(x))

    cdf <- mat.or.vec(r,1)
    cdf[1:r] <- stdn_cdf((x[1:r] - m[1:r])/sqrt(v[1:r]))
    return(cdf)
}
