# Phi computes the standard normal distribution function value at x

Phi <- function(x) {
    y <- 0.5*(1+erf(x/sqrt(2)))
    return(y)
}
