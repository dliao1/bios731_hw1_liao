library(tidyverse)
library(doParallel)
library(foreach)



gen_data <- function(n, beta_true, err_type) {
  beta0 <- 1
  beta_treat <- beta_true
  x <- rbinom(n, 1, prob = 0.5)
  
  while (length(unique(x)) == 1) {
    x <- rbinom(n, 1, prob = 0.5)
  }
  epsilon <- rep(NA, n)
  
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

fit_model <- function(data) {
  model <- lm(y ~ x, data = data)
  return (model)
}

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

# We are given a df of x and y with n rows (for each iteration of the 475 n_sims)
extract_estim_boot_percent <- function(all_boot_betas, beta_true, alpha) {
   # Compute summary statistics
  mean_beta_hat <- mean(all_boot_betas, na.rm = TRUE)  # Mean of bootstrap estimates
  # Percentile confidence interval
  percentile_ci <- quantile(all_boot_betas, probs = c(alpha/2, 1 - (alpha/2)), na.rm = TRUE)
  se_beta_hat = sd(all_boot_betas)
  
  boot_percent_df <- tibble(
    mean_beta_hat = mean_beta_hat,
    se_beta_hat = se_beta_hat,
    ci_l = percentile_ci[1],
    ci_u = percentile_ci[2]
  )
  
  
  boot_percent_df <- boot_percent_df %>%
    mutate(coverage = ifelse(!is.na(ci_l) & !is.na(ci_u) & beta_true >= ci_l & beta_true <= ci_u, 1, 0))

  return(boot_percent_df)
  
}

get_boot_data <- function (original_data, beta_true, sample_size, alpha, nboot) {
  all_boot_betas <- rep(NA, nboot)
  for (i in 1:nboot) {
    boot_sample <- slice_sample(original_data, n = sample_size, replace = TRUE)
    boot_sample_model <- fit_model(data = boot_sample)
    # Extract beta_treatment
    model_estims <- extract_estims(boot_sample_model, beta_true, alpha)
    all_boot_betas[i] <- ifelse(nrow(model_estims) == 1, model_estims$beta_hat, NA)
  }
  return (all_boot_betas)
}


get_boot_t_data <- function(original_data, beta_true, sample_size, alpha, nboot, nboot_t) {
  all_boot_betas <- rep(NA, nboot)
  all_nested_boot_betas <- rep(NA, nboot_t)
  se_stars <- rep(NA, nboot)
  
  
  for (i in 1:nboot) {
    boot_sample <- slice_sample(original_data, n = sample_size, replace = TRUE)
    boot_sample_model <- fit_model(data = boot_sample)
    
    # Extract beta_treatment
    first_model_estims <- extract_estims(boot_sample_model, beta_true, alpha)
    all_boot_betas[i] <- ifelse(nrow(first_model_estims) == 1, first_model_estims$beta_hat, NA)
    
    for (j in 1:nboot_t) {
      nested_boot_sample <- slice_sample(boot_sample, n = sample_size, replace = TRUE)
      nested_boot_sample_model <- fit_model(data = nested_boot_sample)
      nested_model_estims <- extract_estims(model = nested_boot_sample_model, 
                                            beta_true = beta_true, 
                                            alpha = alpha)
      
      all_nested_boot_betas[j] <- ifelse(nrow(nested_model_estims) == 1, nested_model_estims$beta_hat, NA)
      
    }
    se_stars[i] = sd(all_nested_boot_betas,  na.rm = TRUE) 
  }
  
  return (list(boot_betas = all_boot_betas, se_stars = se_stars))
}

extract_estim_boot_t <- function(original_data, all_boot_betas, se_stars, beta_true, alpha) {
  #Original model estimates
  
  original_data_model <- fit_model(data = original_data)
  original_model_estims <- extract_estims(original_data_model, beta_true, alpha)
  beta_hat <- ifelse(nrow(original_model_estims) == 1, original_model_estims$beta_hat, NA)
  
  
  # Compute summary statistics
  mean_beta_hat <- mean(all_boot_betas, na.rm = TRUE)  # Mean of bootstrap estimates
  
  # Percentile confidence interval
  percentile_ci <- quantile(all_boot_betas, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  
  # quantile CI for bootstrap t
  t_stars <- (all_boot_betas - beta_hat) / se_stars
  t_quants = quantile(t_stars, probs = c(alpha/2, 1-(alpha/2)),  na.rm = TRUE)
  se_beta_hat = sd(all_boot_betas)
  
  # lower CI
  boot_t_ci_l <- beta_hat - t_quants[2] * se_beta_hat
  
  # upper CI
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

