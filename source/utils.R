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
  
  boot_percent_df <- tibble(
    mean_beta_hat = mean_beta_hat,
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
    # Non-parametric bootstrap sample
    rows <- sample(1:sample_size, sample_size, replace = TRUE)
    boot_sample <- original_data[rows,]
    boot_sample_model <- fit_model(data = boot_sample)
    # Extract beta_treatment
    model_estims <- extract_estims(boot_sample_model, beta_true, alpha)
    all_boot_betas[i] <- ifelse(nrow(model_estims) == 1, model_estims$beta_hat, NA)
  }
  return (all_boot_betas)
}
  
  
# We are given a df of x and y with n rows (for each iteration of the 475 n_sims)
extract_estim_boot_t <- function(data, beta_true, alpha) {
  original_model_estims <- extract_estims(data, beta_true, alpha)
  beta_hat <- original_model_estims$beta_hat
  
  nboot <- 50  # Number of bootstrap samples
  nboot_t <- 10 # Inner loop to calculate t star
  
  all_boot_betas <- rep(NA, nboot)  # change to rep NA first
  all_nested_boot_betas <- rep(NA, nboot_t) 
  t_star <- rep(NA, nboot_t) 
  
  for (i in 1:nboot) {
    # Non-parametric bootstrap sample
    boot_sample <- slice_sample(data, n = nrow(data), replace = TRUE)

    # Fit linear model and extract beta_treatment
    first_model_estims <- extract_estims(data = boot_sample, 
                                   beta_true = beta_true, 
                                   alpha = alpha)
    
    # Return the estimate or NA if invalid
    all_boot_betas[i] <- ifelse(nrow(first_model_estims) == 1, first_model_estims$beta_hat, NA)
    
    for (j in 1:nboot_t) {
      nested_boot_sample <- slice_sample(boot_sample, n = nrow(boot_sample), replace = TRUE)
      nested_model_estims <- extract_estims(data = nested_boot_sample, 
                                     beta_true = beta_true, 
                                     alpha = alpha)
      
      all_nested_boot_betas[j] <- ifelse(nrow(nested_model_estims) == 1, nested_model_estims$beta_hat, NA)
      
    }
    se_star = sd(all_nested_boot_betas,  na.rm = TRUE) # is this right... idk
    
    if (any(is.na(c(all_boot_betas[i], beta_hat, se_star)))) {
      t_star[i] <- NA
    } else {
      t_star[i] <- (all_boot_betas[i] - beta_hat) / se_star
    }
  }
  
  # Compute summary statistics
  mean_beta_hat <- mean(all_boot_betas, na.rm = TRUE)  # Mean of bootstrap estimates

  # Percentile confidence interval
  percentile_ci <- quantile(all_boot_betas, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  
  # quantile CI for bootstrap t
  t_quants = quantile(t_star, probs = c(alpha/2, 1-(alpha/2)),  na.rm = TRUE)
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
    t_ci_l = boot_t_ci_l, #calculated froim ses of the nested bootstrap
    t_ci_u = boot_t_ci_u,
  )
  
  boot_t_df <- boot_t_df %>%
    mutate(coverage = ifelse(!is.na(percent_ci_l) & !is.na(percent_ci_u) & beta_true >= percent_ci_l & beta_true <= percent_ci_u, 1, 0)) %>%
    mutate(t_coverage = ifelse(!is.na(t_ci_l) & !is.na(t_ci_u) & beta_true >= t_ci_l & beta_true <= t_ci_u, 1, 0))
  return(boot_t_df)
}
