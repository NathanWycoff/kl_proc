
#'A numerically stable softmax function
softmax <- function(x) {
    exp(x - max(x) - log(sum(exp(x - max(x)))))
}

#' (log) L2 distance between two (univariate) Gaussian Mixtures
gmm_l2 <- function(PI1, PI2, mu_1, mu_2, sigma_1, sigma_2) {
    K1 <- length(PI1)
    K2 <- length(PI2)

    ldist <- 0

    # First dist integrated against itself
    for (ki in 1:K1) {
        for (kj in 1:K1) {
            sig <- sqrt(sigma_1[ki]^2+sigma_1[kj]^2)
            ldist <- ldist + PI1[ki] * PI1[kj] * 
                dnorm(mu_1[ki], mu_1[kj], sig)
        }
    }
    # Same with second
    for (ki in 1:K2) {
        for (kj in 1:K2) {
            sig <- sqrt(sigma_2[ki]^2+sigma_2[kj]^2)
            ldist <- ldist + PI2[ki] * PI2[kj] * 
                dnorm(mu_2[ki], mu_2[kj], sig)
        }
    }
    # The cross
    for (ki in 1:K1) {
        for (kj in 1:K2) {
            sig <- sqrt(sigma_1[ki]^2+sigma_2[kj]^2)
            ldist <- ldist - 2 * PI1[ki] * PI2[kj] * 
                dnorm(mu_1[ki], mu_2[kj], sig)
        }
    }

    return(ldist)
}

rdirich <- function(alpha) {
    x <- rgamma(length(alpha), alpha, 1)
    x / sum(x)
}

#' Variance of a univariate normal mixture model
#' @param mu The vector of means
#' @param sigma The vector of standard deviations
#' @param PI The vector of mixture weights
#' @return The scalar variance.
gmm_var <- function(mu, sigma, PI) {
    N <- length(mu)
    var_term <- crossprod(PI, sigma^2)
    Z_VAR <- matrix(NA, nrow = N, ncol = N)
    for (i in 1:N) {
        for (j in 1:i) {
            if (i == j) {
                Z_VAR[i,i] <- PI[i] * (1-PI[i])
            } else {
                Z_VAR[i,j] <- Z_VAR[j,i] <- -PI[i] * PI[j]
            }
        }
    }

    mean_term <- t(mu) %*% Z_VAR %*% mu
    return(as.numeric(var_term + mean_term))
}
