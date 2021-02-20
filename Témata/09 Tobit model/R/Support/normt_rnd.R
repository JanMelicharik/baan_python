# PURPOSE: random draws from a normal truncated to (left,right) interval
# ------------------------------------------------------
# USAGE: y = normt_rnd(mu,sigma2,left,right)
# where:   mu = mean (nobs x 1)
#      sigma2 = variance (nobs x 1)
#        left = left truncation points (nobs x 1)
#       right = right truncation points (nobs x 1)
# ------------------------------------------------------
# RETURNS: y = (nobs x 1) vector
# ------------------------------------------------------
# NOTES: use y = normt_rnd(mu,sigma2,left,mu+5*sigma2)
#        to produce a left-truncated draw
#        use y = normt_rnd(mu,sigma2,mu-5*sigma2,right)
#        to produce a right-truncated draw
# ------------------------------------------------------
# SEE ALSO: normlt_rnd (left truncated draws), normrt_rnd (right truncated)
#
#
# adopted from Bayes Toolbox by
# James P. LeSage, Dept of Economics
# University of Toledo
# 2801 W. Bancroft St,
# Toledo, OH 43606
# jpl@jpl.econ.utoledo.edu

normt_rnd <- function(mu,sigma2,left,right) {
    std <- sqrt(sigma2)

    points_left <- which(left == -999)
    points_right <- which(right == 999)
    a_term <- norm_cdf((left - mu)/std)
    a_term[points_left] <- 0

    b_term <- norm_cdf((right - mu)/std)
    b_term[points_right] <- 1

    uniforms <- runif(length(mu))

    p <- a_term + as.vector(uniforms*(b_term - a_term))
    # Calculate bounds on probabilities
    #lowerProb <- Phi((left-mu)/std)
    #upperProb <- Phi((right-mu)/std)
    # Draw uniform from within (lowerProb,upperProb)
    #u <- lowerProb + (upperProb - lowerProb)*runif(length(mu))
    # Find needed quantiles
    result <- mu + std*norm_inv(as.matrix(p))
    return(result)
}
