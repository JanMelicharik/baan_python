# PURPOSE: returns the pdf at x of the gamma(a) distribution
# ---------------------------------------------------
# USAGE: pdf = gamm_pdf(x,a)
#  where: x = a vector
#         a = a scalar for gamma(a)
#---------------------------------------------------
#  RETURNS:
#    a vector of pdf at each element of x of the gamma(a) distribution
# --------------------------------------------------
#
#     Anders Holtsberg, 18-11-93
#     Copyright (c) Anders Holtsberg
gamm_pdf <- function(x,a) {
    f <- x^(a - 1)*exp(-x)/gamma(a)
    I0 <- which(x < 0)
    f[I0] <- 0
    return(f)
}
