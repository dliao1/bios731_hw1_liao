library(tidyverse)
library(doParallel)
library(foreach)



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

# We are given a df of x and y with n rows (for each iteration of the 475 n_sims)
extract_estim_boot <- function(simdata, beta_true, alpha) {
  nboot <- 1  # Number of bootstrap samples, 1 right now so its fast
  all_boot_estims <- numeric(nboot)  # Preallocate a numeric vector
  
  for (b in 1:nboot) {
    # Non-parametric bootstrap sample
    indices <- sample(seq_len(nrow(simdata)), size = nrow(simdata), replace = TRUE)
    boot_sample <- simdata[indices, , drop = FALSE]
    
    # Fit linear model
    model <- lm(y ~ x, data = boot_sample)
    
    # Extract beta_treatment
    beta_hat_wrapper <- tidy(model, conf.int = TRUE) %>%
      filter(term == "x") %>%
      rename(beta_hat = estimate) %>%
      select(beta_hat)
    
    # Return the estimate or NA if invalid
    if (nrow(beta_hat_wrapper) == 1) {
      all_boot_estims[b] <- beta_hat_wrapper$beta_hat
    } else {
      all_boot_estims[b] <- NA
    }
  }
  
  # Compute summary statistics
  mean_beta_hat <- mean(all_boot_estims, na.rm = TRUE)  # Mean of bootstrap estimates
  se_beta_hat <- sd(all_boot_estims, na.rm = TRUE) / sqrt(length(all_boot_estims))  # Standard error?? (this is for Wald type so not what we want but its a start)
  
  # Percentile confidence interval
  percentile_ci <- quantile(all_boot_estims, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  
  tibble(
    mean_beta_hat = mean_beta_hat,
    se_beta_hat = se_beta_hat,
    ci_l = percentile_ci[1],
    ci_u = percentile_ci[2]
  )
}
