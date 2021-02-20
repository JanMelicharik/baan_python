#  PURPOSE: returns the cdf at x of the gamma(a) distribution
# --------------------------------------------------
#  USAGE: cdf = gamm_cdf(x,a)
#  where: x = a vector
#         a = a scalar gamma(a)
# ---------------------------------------------------
#  RETURNS:
#     a vector of cdf at each element of x of the gamma(a) distribution
#  --------------------------------------------------
#
#        Anders Holtsberg, 18-11-93
#        Copyright (c) Anders Holtsberg

gamm_cdf <- function(x,a) {
    b <- unname(gammainc(x,a))
    cdf <- b[3]
    I0 <- which(x < 0)
    cdf[I0] <- 0
    return(cdf)
}
