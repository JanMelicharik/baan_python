# PURPOSE: computes the standard normal probability density
#          for each component of x with mean=0, variance=1
# ---------------------------------------------------
# USAGE: cdf = stdn_pdf(x)
# where: x = (n x m) matrix
# ---------------------------------------------------
# RETURNS: pdf = (nxm) matrix containing pdf
#                  for each element of x
# ---------------------------------------------------
#
# Written by TT (Teresa.Twaroch@ci.tuwien.ac.at)
# Updated by KH (Kurt.Hornik@ci.tuwien.ac.at)

stdn_pdf <- function(x) {
    r <- nrow(x)
    c <- ncol(x)
    s <- r*c
    x <- matrix(x,1,s)
    pdf <- mat.or.vec(1,s)

    k <- which(is.na(x))
    if(any(k)) {
        pdf[k] <- NA
    }

    k <- which(!is.infinite(x))
    if(any(k)) {
        pdf[k] <- (2*pi)^(-1/2)*exp(-x[k]^2/2)
    }
    pdf <- matrix(pdf,r,c)
    return(pdf)
}
