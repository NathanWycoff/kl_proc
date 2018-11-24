#!/usr/bin/Rscript
#  R/lin_reg.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 11.23.2018

## Apply our stochastic process to linear data.
source('R/stoch_process.R')
set.seed(123)
N <- 10
x <- seq(0, 1, length.out = N)

y <- 2*x + rnorm(N,0,0.1)

eps <- 2e-3
DELTA <- eps * as.matrix(dist(x))


plot(x, y)
