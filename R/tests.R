#!/usr/bin/Rscript
#  R/tests.R Author "Nathan Wycoff <nathanbrwycoff@gmail.com>" Date 11.23.2018

source("R/lib.R")

## Test that L2 distance works
# Distance to yourself is 0 
gmm_l2(PI1 = c(0.5,0.5), PI2 = c(0.5,0.5), mu_1 = c(10,1), mu_2 = c(10,1),
                   sigma_1 = c(20,0.1), sigma_2 = c(20, 0.1))
gmm_l2(PI1 = c(0.5,0.5), PI2 = c(0.5,0.5), mu_1 = c(1,10), mu_2 = c(10,1),
                   sigma_1 = c(20,0.1), sigma_2 = c(0.1, 20))

# Distance to someone else is not 0
gmm_l2(PI1 = c(0.5,0.5), PI2 = c(0.1,0.9), mu_1 = c(1,10), mu_2 = c(10,1),
                   sigma_1 = c(20,0.1), sigma_2 = c(0.1, 20))
