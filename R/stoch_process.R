#!/usr/bin/Rscript
#  R/stoch_process.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 11.23.2018

## The entire enchilada, the Bound Stochastic Process
get_count <- function(zs, K) {
    counts <- rep(0, K)
    for (l in 1:length(zs)) {
        counts[zs[l]] <- counts[zs[l]] + 1
    }
    return(counts)
}

#' Do GMM inference, constraining the estimate to be within some distance of the truth
#' @param Xs A list of observed data, each element of the list gives observations at a location. Some list elements may have zero observations (corresponds to locations where predictions are desired).
#' @param K The number of mixture components
#' @param Zs The mixture assignments corresponding to each of the Xs.
#' @param mu_init The initializations for the means.
#' @param sigmas The fixed standard deviations (length K).
#' @param Delta The length(Xs) by length(Xs) matrix giving maximum allowed pairwise deviation in terms of L2 norm in estimated distribution at any point in the MCMC.
bound_gmm <- function(K, Xs, Zs, mu_init, sigmas, DELTA, xx) {
    G <- length(Xs)
    for (iter in 1:iters) {
        for (g in 1:G) {
            xs <- Xs[[g]]
            if (length(xs) > 0) {

        }
    }
}
