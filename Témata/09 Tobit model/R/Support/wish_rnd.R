wish_rnd <- function(sigma,v) {
    n <- nrow(sigma)
    k <- ncol(sigma)

    if(n != k) {
        stop('wish_rnd: requires a square matrix')
    } else if(n < k) {
        stop('wish_rnd: n must be >= k+1 for a finite distribution')
    }
    t <- chol(sigma)
    p <- eigen(t)
    if(sum(p$values < 0) > 0) {
        stop('wish_rnd: matrix must be a positive definite')
    }
    a <- matrix(0,n,v)
    for(ii in 1:n) {
        a[ii,] <- rnorm(v)
    }
    y <- t(t) %*% a #matrix(rnorm(),n,v)
    w <- y  %*% t(y)
    return(w)
}
