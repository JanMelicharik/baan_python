# PURPOSE: compute random draws from a left-truncated normal
#           distribution, with mean = mu, variance = sigma2
#  ------------------------------------------------------
# USAGE: y = normlt_rnd(mu,sigma2,left)
# where:   mu = mean (scalar or vector)
#       sigma2 = variance (scalar or vector)
#         left = left truncation point (scalar or vector)
# ------------------------------------------------------
# RETURNS: y = (scalar or vector) the size of mu, sigma2
# ------------------------------------------------------
# NOTES: This is merely a convenience function that
#       calls normt_rnd with the appropriate arguments
# ------------------------------------------------------
#
# written by:
# James P. LeSage, Dept of Economics
# University of Toledo
# 2801 W. Bancroft St,
# Toledo, OH 43606
# jpl@jpl.econ.utoledo.edu

normlt_rnd <- function(mu,sigma2,left) {
    nobs <- length(mu)
    right <- 999*rep(1,nobs)

    result <- normt_rnd(mu,sigma2,left,right)
    return(result)
}
