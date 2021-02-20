# PURPOSE: random multivariate random vector based on
#          var-cov matrix sig
#---------------------------------------------------
# USAGE:   y = norm_rnd(sig)
# where:   sig = a square-symmetric covariance matrix
# NOTE: for mean b, var-cov sig use: b +  norm_rnd(sig)
#---------------------------------------------------
# RETURNS: y = random vector normal draw mean 0, var-cov(sig)
 #---------------------------------------------------

# written by:
# James P. LeSage, Dept of Economics
# University of Toledo
# 2801 W. Bancroft St,
# Toledo, OH 43606
# jlesage@spatial-econometrics.com

norm_rnd <- function(sig){
    h = chol(sig)
    rv = rnorm(nrow(as.matrix(sig)))
    res = t(h) %*% rv

    return(res)
}
