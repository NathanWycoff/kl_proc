#!/usr/bin/Rscript
#  R/playground.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 11.23.2018

source('R/lib.R')

#' Sample from a gaussian mixture model.
rgmm <- function(N, K, mus, sigmas, PI) {
    xs <- rep(NA, N)
    zs <- rep(NA, N)
    for (n in 1:N) {
        zs[n] <- sample(1:K, 1, prob = PI)
        xs[n] <- rnorm(1,mus[zs[n]], sigmas[zs[n]])
    }
    ret <- list(x = xs, z = zs)
}

get_count <- function(zs, K) {
    counts <- rep(0, K)
    for (l in 1:length(zs)) {
        counts[zs[l]] <- counts[zs[l]] + 1
    }
    return(counts)
}

# Sample, with the constraint that we don't go eps from the true distribution.
#eps <- 1e-3
get_dist <- function(PI_hat, mu_hat, sigma_hat) {
    gmm_l2(PI1 = PI_hat, PI2 = PI, mu_1 = mu_hat, mu_2 = mus,
                       sigma_1 = sigmas, sigma_2 = sigma_hat)
}

# Do GMM inference, constraining the estimate to be within some distance of the truth
bound_gmm <- function(x, z_init, mu_init, sigmas, eps, xx) {
    z_hat <- z_init
    mu_hat <- mu_init
    pi_hat <- (get_count(z_hat, K) + alpha) / (length(z_hat) + K*alpha)

    for (iter in 1:iters) {
        print(iter)

        # Iterate through each point and pick a new cluster assignment
        for (n in 1:N) {
            # All points except this one.
            x_me <- ret$x[n]
            x_other <- ret$x[-n]

            # Remove current point from cluster proportion calculations
            z_hat_arch <- z_hat
            z_hat_arch[n] <- -1
            z_hat <- z_hat[-n]

            # Sample new cluster asignments 
            logprobs <- sapply(1:K, function(k) {
                        members <- log(pi_hat[k])
                        density <- dnorm(x_me, mu_hat[k], sigmas[k], log = TRUE)
                        return(members + density)
                })
            probs <- softmax(logprobs)

            # Sample new point
            asgn <- sample(1:K, 1, prob = probs)
            z_hat_arch[n] <- asgn
            z_hat <- z_hat_arch
        }

        # Update Mixture Weights
        clust_counts <- get_count(z_hat, K) + alpha
        dist <- Inf
        while (dist > eps) {
            pi_hat <- rdirich(clust_counts)
            dist <- get_dist(pi_hat, mu_hat, sigmas)
        }

        # Sample new params
        for (k in 1:K) { 
            post_var <- 1/(1/mu_sigma^2 + sum(z_hat==k) / sigmas[k]^2)
            post_mean <- post_var * (mu_mu / mu_sigma^2 + sum(ret$x[z_hat==k]) / sigmas[k]^2)
            mu_maybe <- mu_hat
            dist <- Inf
            while (dist > eps) {
                print("Oh we did a thing!")
                proposal <- rnorm(1, post_mean, sqrt(post_var))
                mu_maybe[k] <- proposal
                dist <- get_dist(pi_hat, mu_maybe, sigmas)
            }
            mu_hat <- mu_maybe
        }

        mu_trace[iter,] <- mu_hat
        PI_trace[iter,] <- (get_count(z_hat, K) + alpha) / (length(z_hat) + K*alpha)
    }

    dense_trace <- matrix(NA, nrow = iters, ncol = length(xx))
    var_trace <- matrix(NA, nrow = iters, ncol = length(xx))
    for (ix in 1:length(xx)) {
        for (iter in 1:iters) {
            dens_vec <- sapply(1:K, function(k) dnorm(xx[ix], mu_trace[iter,k], sigmas[k]))
            dense_trace[iter, ix] <- dens_vec %*% PI_trace[iter,]
            var_trace[iter, ix] <- gmm_var(mu_trace[iter,], sigmas, PI_trace[iter,])
        }
    }

    pred_dens <- colMeans(dense_trace)
    pred_var <- colMeans(var_trace)

    return(list(pred_dens, pred_var))
}
