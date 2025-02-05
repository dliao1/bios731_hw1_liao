library(tidyverse)
library(doParallel)
library(foreach)
library(tictoc)

# Extract estimates from linear model
extract_estims <- function(model, beta_true, alpha) {
  estims_df <- tidy(model, conf.int = TRUE) %>%
    filter(term == "x") %>%
    rename(beta_hat = estimate) %>%
    rename(se_beta = std.error) %>%
    mutate(beta_diff = beta_hat - beta_true) %>%
    mutate(coverage = ifelse(beta_true >= conf.low & beta_true <= conf.high, 1, 0)) %>%
    select(beta_hat, beta_diff, se_beta, coverage)
  return (estims_df)
}

# Extracts estimates and bootstrap percentile intervals
extract_estim_boot_percent <- function(all_boot_betas, beta_true, alpha) {
  # Mean of bootstrap estimates
  mean_beta_hat <- mean(all_boot_betas, na.rm = TRUE)  
  
  # Percentile confidence interval
  percentile_ci <- quantile(all_boot_betas, probs = c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)
  
  # Calculates se_beta_hat
  se_beta_hat = sd(all_boot_betas)
  
  boot_percent_df <- tibble(
    mean_beta_hat = mean_beta_hat,
    se_beta_hat = se_beta_hat,
    ci_l = percentile_ci[1],
    ci_u = percentile_ci[2]
  )
  
  # Calculates coverage
  boot_percent_df <- boot_percent_df %>%
    mutate(coverage = ifelse(!is.na(ci_l) & !is.na(ci_u) & beta_true >= ci_l & beta_true <= ci_u, 1, 0))
  
  return(boot_percent_df)
  
}

# Extracts estimates and bootstrap t intervals
extract_estim_boot_t <- function(original_data, all_boot_betas, se_stars, beta_true, alpha) {
  
  #Original model estimates
  original_data_model <- fit_model(data = original_data)
  original_model_estims <- extract_estims(original_data_model, beta_true, alpha)
  beta_hat <- ifelse(nrow(original_model_estims) == 1, original_model_estims$beta_hat, NA)
  
  # Mean of bootstrap estimates
  mean_beta_hat <- mean(all_boot_betas, na.rm = TRUE)
  
  # Percentile confidence interval
  percentile_ci <- quantile(all_boot_betas, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  
  # Bootstrap t interval
  t_stars <- (all_boot_betas - beta_hat) / se_stars
  t_quants = quantile(t_stars, probs = c(alpha/2, 1-(alpha/2)),  na.rm = TRUE)
  
  # Calculates se_beta_hat
  se_beta_hat = sd(all_boot_betas)
  
  # Lower CI
  boot_t_ci_l <- beta_hat - t_quants[2] * se_beta_hat # se_beta_hat estimated from top level bootstrap
  
  # Upper CI
  boot_t_ci_u <- beta_hat - t_quants[1] * se_beta_hat
  
  
  boot_t_df <- tibble(
    mean_beta_hat = mean_beta_hat,
    se_beta_hat = se_beta_hat, # same as the percentile t method
    percent_ci_l = percentile_ci[1],
    percent_ci_u = percentile_ci[2],
    t_ci_l = boot_t_ci_l, #calculated from ses of the nested bootstrap
    t_ci_u = boot_t_ci_u,
  )
  
  boot_t_df <- boot_t_df %>%
    mutate(coverage = ifelse(!is.na(t_ci_l) & !is.na(t_ci_u) & beta_true >= t_ci_l & beta_true <= t_ci_u, 1, 0))
  return(boot_t_df)
  
}