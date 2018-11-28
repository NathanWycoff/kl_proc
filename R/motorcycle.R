#!/usr/bin/Rscript
#  R/motorcycle.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 11.28.2018

## Can our model capture the heteroscedasticity in the famous motorcycle data?
require(MASS)
source('R/stoch_process.R')

df <- mcycle

# Get rid of rows with less than some distance, as we can't handle small/zero distance yet.
d <- as.matrix(dist(df$times))
diag(d) <- Inf
min_dist <- 0.01
df <- df[apply(d, 1, min) > min_dist,]

# Standardize
df$accel <- (df$accel - mean(df$accel)) / sd(df$accel)

# Get variance
eps <- 0.2
DELTA <- eps * as.matrix(dist(df$times))

K <- 3
Ys <- as.list(df$accel)
Zs <- as.list(sample(1:K, nrow(df), replace = TRUE))
mu_init <- c(-1, 0, 1)
sigmas <- c(0.5,0.5,0.5)

ret <- bound_gmm(K, Ys, Zs, mu_init, sigmas, DELTA, iters = 10, mu_mu = 0)

#plot(df, ylim = c(-100, 100))
plot(df)
lines(df$times, ret$mu, col = 'red')
lines(df$times, ret$mu + 2*ret$sig, col = 'blue')
lines(df$times, ret$mu - 2*ret$sig, col = 'blue')
