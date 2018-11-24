#!/usr/bin/Rscript
#  R/normal_mixture.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 11.23.2018

source('R/playground.R')
## Do Gibbs sampling, for now with variances known/fixed
set.seed(123)

## A normal mixture model, estimated via Gibbs sampling.
K <- 3
mus <- c(-4,2,8)
sigmas <- c(3,3,3)
PI <- c(1/3,1/3,1/3)
N <- 1

# Predictive locations
xx <- seq(-10,10, length.out = 1e3)

# Some hyperparams
mu_mu <- 0
mu_sigma <- 10

ret <- rgmm(N, K, mus, sigmas, PI)

# Smoothing param so that we can visit empty clusters, equiv to exchangible dirichlet with prob alpha
alpha <- 1

iters <- 100
sigma <- 1

# Estimated cluster assignemnts
# init at truth for now
z_hat <- sample(1:K, N, replace = TRUE)
mu_hat <- mus

mu_trace <- matrix(NA, nrow = iters, ncol = K)
PI_trace <- matrix(NA, nrow = iters, ncol = K)

set.seed(123)

xx <- seq(-10, 10, length.out = 100)
eps <- 1.0e-4
pred_dens <- bound_gmm(ret$x, ret$z, mus, sigmas, eps, xx)

hist(ret$x, probability = TRUE, xlim = c(-10, 10), breaks = 25, col ='skyblue')
points(xx, pred_dens, type = 'l', col = 'blue', lwd = 5)
