library(tidyverse)
library(doParallel)
library(foreach)
library(tictoc)

# Generates data given n, true beta, and error type (1 for N(0,2) errors, 0 for lognormal(0,2) errors)
gen_data <- function(n, beta_true, err_type) {
  beta0 <- 1
  beta_treat <- beta_true
  x <- rbinom(n, 1, prob = 0.5)
  
  # Makes sure all xs are not the same so linear model fit later won't fail
  while (length(unique(x)) == 1) {
    x <- rbinom(n, 1, prob = 0.5)
  }
  
  epsilon <- rep(NA, n)
  
  # Determines what kind of errors to generate
  if (err_type == 1) {
    epsilon <- rnorm(n, mean = 0, sd = sqrt(2))
  } else {
    epsilon <- rlnorm(n, mean = 0, sd = sqrt(log(2)))
  }
  
  y = beta0 + beta_treat * x + epsilon
  tibble(
    x = x,
    y = y,
  )
}





