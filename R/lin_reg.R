#!/usr/bin/Rscript
#  R/lin_reg.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 11.23.2018

## Apply our stochastic process to linear data.
source('R/stoch_process.R')
N <- 10
x <- seq(0, 1, length.out = N)

y <- 2*x + rnorm(N,0,0.1)

eps <- 16
DELTA <- eps * as.matrix(dist(x))

K <- 3
Ys <- as.list(y)
Zs <- as.list(sample(1:K, N, replace = TRUE))

mu_init <- c(mean(y) - 1, mean(y), mean(y) + 1)
sigmas <- rep(0.1, K)

#set.seed(123)

ret <- bound_gmm(K, Ys, Zs, mu_init, sigmas, DELTA, iters = 10, mu_mu = 0)

plot(x, y, ylim = c(-100, 100))
lines(x, ret$mu, col = 'red')
lines(x, ret$mu + 2*ret$sig, col = 'blue')
lines(x, ret$mu - 2*ret$sig, col = 'blue')
