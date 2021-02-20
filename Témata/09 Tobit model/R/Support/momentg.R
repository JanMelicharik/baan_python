# PURPOSE: computes Gewke's convergence diagnostics NSE and RNE
#          (numerical std error and relative numerical efficiencies)
# -----------------------------------------------------------------
# USAGE: result = momentg(draws)
# where: draws = a matrix of Gibbs draws (ndraws x nvars)
# ----------------------------------------------------------------
#  RETURNS: a structure result:
#           result.meth     = 'momentg'
#           result.ndraw    = # of draws
#           result.nvar     = # of variables
#           result(i).pmean = posterior mean for variable i
#           result(i).pstd  = posterior std deviation
#           result(i).nse   = nse assuming no serial correlation for variable i
#          result(i).rne   = rne assuming no serial correlation for variable i
#           result(i).nse1  = nse using 4% autocovariance tapered estimate
#           result(i).rne1  = rne using 4% autocovariance taper
#           result(i).nse2  = nse using 8% autocovariance taper
#           result(i).rne2  = rne using 8% autocovariance taper
#           result(i).nse3  = nse using 15% autocovariance taper
#           result(i).rne3  = rne using 15% autocovariance taper
#  -----------------------------------------------------------------
#  SEE ALSO: coda(), apm()
#  -----------------------------------------------------------------
#  REFERENCES: Geweke (1992), `Evaluating the accuracy of sampling-based
#  approaches to the calculation of posterior moments', in J.O. Berger,
#  J.M. Bernardo, A.P. Dawid, and A.F.M. Smith (eds.) Proceedings of
#  the Fourth Valencia International Meeting on Bayesian Statistics,
#  pp. 169-194, Oxford University Press
#  Also: `Using simulation methods for Bayesian econometric models:
#  Inference, development and communication', at: www.econ.umn.edu/~bacc
#  -----------------------------------------------------------------
# # % written by:
#  James P. LeSage, Dept of Economics
#  University of Toledo
#  2801 W. Bancroft St,
#  Toledo, OH 43606
#  jpl@jpl.econ.utoledo.edu
#
#  NOTE: this code draws heavily on MATLAB programs written by
#  Siddartha Chib available at: www.econ.umn.edu/~bacc
#  I have repackaged it to make it easier to use.

momentg <- function(draws) {
    ndraw <- nrow(draws)
    nvar <- ncol(draws)

    results <- list()
    results$ndraw <- ndraw
    results$meth <- 'momentg'
    results$nvar <- nvar

    NG <- 100

    if(ndraw < NG) {
        stop('momentg: needs a larger number of ndraws')
    }

    ntaper <- c(4,8,15)
    ns <- floor(ndraw/NG)
    nuse <- ns*NG

    for(jf in 1:nvar) { # loop over all variables
        cnt <- 0
        cn <- mat.or.vec(NG,NG)
        cd <- mat.or.vec(NG,NG)
        cdn <- mat.or.vec(NG,NG)
        cdd <- mat.or.vec(NG,NG)
        cnn <- mat.or.vec(NG,NG)
        cvar <- mat.or.vec(NG,NG)

        # form sufficiency statistics needed below
        td <- 0
        tn <- 0
        tdd <- 0
        tnn <- 0
        tdn <- 0
        tvar <- 0

        for(ig in 1:NG) {
            gd <- 0
            gn <- 0
            gdd <- 0
            gdn <- 0
            gnn <- 0
            gvar <- 0

            for(is in 1:ns) {
                cnt <- cnt + 1
                g <- draws[cnt,jf]
                ad <- 1
                an <- ad*g
                gd <- gd + ad
                gn <- gn + an
                gdn <- gdn + ad*an
                gdd <- gdd + ad*ad
                gnn <- gnn + an*an
                gvar <- gvar + an*g
            }

            td <- td + gd
            tn <- tn + gn
            tdn <- tdn + gdn
            tdd <- tdd + gdd
            tnn <- tnn + gnn
            tvar <- tvar + gvar

            cn[ig] = gn/ns
            cd[ig] = gd/ns
            cdn[ig] = gdn/ns
            cdd[ig] = gdd/ns
            cnn[ig] = gnn/ns
            cvar[ig] = gvar/ns
        }

        eg <- tn/td
        varg <- tvar/td - eg^2
        sdg <- -1

        if(varg > 0) {
            sdg <- sqrt(varg)
        }

        # save posterior means and std deviations to results structure
        results$pmean[jf] <- eg
        results$pstd[jf] <- sdg

        # numerical standard error assuming no serial correlation
        varnum <- (tnn - 2*eg*tdn + tdd*eg^2)/(td^2)
        sdnum <- -1

        if(varnum > 0) {
            sdnum <- sqrt(varnum)
        }

        # save to results structure
        results$nse[jf] <- sdnum
        results$rne[jf] <- varg/(nuse*varnum)

        # get autocovariance of grouped means
        barn <- tn/nuse
        bard <- td/nuse

        for(ig in 1:NG) {
            cn[ig] <- cn[ig] - barn
            cd[ig] <- cd[ig] - bard
        }

        rnn <- 0
        rdd <- 0
        rnd <- 0
        rdn <- 0

        for(lag in 0:(NG-1)) {
            ann <- 0
            add <- 0
            and <- 0
            adn <- 0
            for(ig in (lag+1):NG) {
                ann <- ann + cn[ig]*cn[ig-lag]
                add <- add + cd[ig]*cd[ig-lag]
                and <- and + cn[ig]*cd[ig-lag]
                adn <- adn + cd[ig]*cn[ig]
            }
            # index 0 not allowed, lag+1 stands for lag
            rnn[lag+1] <- ann/NG
            rdd[lag+1] <- add/NG
            rnd[lag+1] <- and/NG
            rdn[lag+1] <- adn/NG
        }

        # numerical standard error with tapered autocovariance functions
        for(mm in 1:3) {
            m <- ntaper[mm]
            am <- m
            snn <- rnn[1]
            sdd <- rdd[1]
            snd <- rnd[1]
            for(lag in 1:(m-1)) {
                att <- 1 - lag/am
                snn <- snn + 2*att*rnn[lag+1]
                sdd <- sdd + 2*att*rdd[lag+1]
                snd <- snd + att*(rnd[lag+1] + rnd[lag+1])
            }

            varnum <- ns*nuse*(snn - 2*eg*snd + sdd*eg^2)/(td^2)
            sdnum <- -1

            if(varnum > 0) {
                sdnum <- sqrt(varnum)
            }

            # save results to list
            if(mm == 1) {
                results$nse1[jf] <- sdnum
                results$rne1[jf] <- varg/(nuse*varnum)
            } else if(mm == 2) {
                results$nse2[jf] <- sdnum
                results$rne2[jf] <- varg/(nuse*varnum)
            } else if(mm == 3) {
                results$nse3[jf] <- sdnum
                results$rne3[jf] <- varg/(nuse*varnum)
            }
        }
    }
    return(results)
}
