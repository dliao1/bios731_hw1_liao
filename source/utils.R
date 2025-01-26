library(tidyverse)


gen_data <- function(n, beta_true, err_type) {
  beta0 <- 1
  beta_treat <- beta_true
  x <- rbinom(n, 1, prob = 0.5)
  epsilon <- ifelse (err_type == 1, rnorm(n, mean = 0, sd = sqrt(2)), rlnorm(n, mean = 0, sd = sqrt(log(2))))
  y = beta0 + beta_treat * x + epsilon
  
  tibble(
    x = x,
    y = y
  )
}

extract_estim_wald <- function(simdata, beta_true, alpha) {
  x <- simdata$x
  y <- simdata$y
  model <- summary(lm(y ~ x))
  
  estims_df <- tidy(model, conf.int = TRUE) %>%
    filter(term == "x") %>%
    mutate(coverage = ifelse(beta_true >= conf.low & beta_true <= conf.high, 1, 0)) %>%
    rename(beta_hat = estimate) %>%
    rename(se_beta = std.error) %>%
    mutate(beta_diff = beta_hat - beta_true) %>%
    mutate(ci_l = beta_hat - qnorm(1 - alpha/ 2) * se_beta) %>%
    mutate(ci_u = beta_hat + qnorm(1 - alpha/ 2) * se_beta) %>%
    select(beta_hat, beta_diff, ci_l, ci_u, se_beta, coverage)
  
 
  return (estims_df)
  
  
}