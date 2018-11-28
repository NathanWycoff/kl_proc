#!/usr/bin/Rscript
#  R/stoch_process.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 11.23.2018
source("R/lib.R")

## The entire enchilada, the Bound Stochastic Process
get_count <- function(zs, K) {
    counts <- rep(0, K)
    for (l in 1:length(zs)) {
        counts[zs[l]] <- counts[zs[l]] + 1
    }
    return(counts)
}

#' Do GMM inference, constraining the estimate to be within some distance of the truth
#' @param Ys A list of observed data, each element of the list gives observations at a location. Some list elements may have zero observations (corresponds to locations where predictions are desired).
#' @param K The number of mixture components
#' @param Zs The mixture assignments corresponding to each of the Ys.
#' @param mu_init The initializations for the means.
#' @param sigmas The fixed standard deviations (length K).
#' @param Delta The length(Ys) by length(Ys) matrix giving maximum allowed pairwise deviation in terms of L2 norm in estimated distribution at any point in the MCMC.
#' @param alpha Exchangible dirichlet prior on mixture weights
#' @param mu_mu The prior mean on the expectation of each normal mixture component
#' @param mu_sigma The prior standard deviation on the expectation of each normal mixture component
#' @return A list with mean and standard deviation at each location.
bound_gmm <- function(K, Ys, Zs, mu_init, sigmas, DELTA, iters = 100, alpha = 1, mu_mu = 0,
                      mu_sigma = 100, verbose = 0) {
    # Default params for debugging purposes (this needs to be commented out unless debugging
    #iters <- 10
    #alpha <- 1
    #mu_mu <- 0
    #mu_sigma <- 10

    ### Initialize parameter matrices
    G <- length(Ys)
    MU <- matrix(mu_init, nrow = G, ncol = K, byrow = TRUE)

    zs <- unlist(Zs)
    pi_init <- (get_count(zs, K) + alpha) / (length(zs) + K*alpha)
    PI <- matrix(pi_init, nrow = G, ncol = K, byrow = TRUE)

    pred_mean_trace <- matrix(NA, nrow = iters, ncol = G)
    pred_var_trace <- matrix(NA, nrow = iters, ncol = G)
    total_samples <- 0# Number of times PI/MU are sampled

    for (iter in 1:iters) {
        if (verbose >= 0) {

            print(iter)
        }
        for (g in 1:G) {
            if (verbose >= 1) {
                print("Starting Point Assignment")
            }
            ys <- Ys[[g]]
            zs <- Zs[[g]]

            # Do cluster assignment for this group
            if (length(ys) > 0) {
                for (n in 1:length(ys)) {
                    if (verbose >= 2) {
                        print(paste("On Data Group", g))
                    }
                    # All points except this one.
                    y_me <- ys[n]
                    y_other <- ys[-n]

                    # Remove current point from cluster proportion calculations
                    zs_arch <- zs
                    zs_arch[n] <- -1
                    zs <- zs[-n]

                    # Sample new cluster asignments 
                    logprobs <- sapply(1:K, function(k)
                                       log(PI[g,k]) + dnorm(y_me, MU[g,k], sigmas[k], log = TRUE))
                    probs <- softmax(logprobs)

                    # Sample new point
                    zs_arch[n] <- sample(1:K, 1, prob = probs)
                    zs <- zs_arch
                }
            }

            # Sample Parameters for Each Mixture, as well as mixture weights
            if (verbose >= 1) {
                print("Starting Parameter Estimation")
            }
            clust_counts <- get_count(zs, K) + alpha
            dists <- rep(Inf, G)
            while (max(dists - DELTA[g,]) > sqrt(.Machine$double.eps)) {
                # Draw mixture weights
                pi_prop <- rdirich(clust_counts)

                # Draw Means
                mu_prop <- rep(NA, K)
                for (k in 1:K) {
                    post_var <- 1/(1/mu_sigma^2 + sum(zs==k) / sigmas[k]^2)
                    post_mean <- post_var * (mu_mu / mu_sigma^2 + sum(ys[zs==k]) / sigmas[k]^2)
                    mu_prop[k] <- rnorm(1, post_mean, sqrt(post_var))
                }

                # Determine distance to other clusters
                dists <- rep(-1, G)
                for (g1 in 1:G) {
                    if (g1 != g) {
                        dists[g1] <- gmm_l2(pi_prop, PI[g1,], mu_prop, MU[g1,], sigmas, sigmas)
                    }
                }

                total_samples <- total_samples + 1
            }

            # Update Params
            MU[g,] <- mu_prop
            PI[g,] <- pi_prop

            # Store predictive mean at Ys
            pred_mean_trace[iter, g] <- t(MU[g,]) %*% PI[g,]
            pred_var_trace[iter, g] <- gmm_var(MU[g,], sigmas, PI[g,])
        }
    }

    print(paste(total_samples / (G*iters), "mean samples per iteration"))

    return(list(mu = colMeans(pred_mean_trace), sig = colMeans(sqrt(pred_var_trace))))
}
