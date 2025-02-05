library(tidyverse)
library(doParallel)
library(foreach)
library(tictoc) 


# Returned value: nboot long list of estimated betas 
get_boot_data <- function (original_data, beta_true, sample_size, alpha, nboot) {
  all_boot_betas <- rep(NA, nboot)
  
  # Draws nboot number of bootstrap samples, fits model from current bootstrap
  # sample, and adds estimated beta to list.
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
  
  # Draws nboot number of bootstrap samples, fits model from current bootstrap
  # sample, and adds estimated beta to list.
  for (i in 1:nboot) {
    boot_sample <- slice_sample(original_data, n = sample_size, replace = TRUE)
    boot_sample_model <- fit_model(data = boot_sample)
    
    # Extract beta_treatment
    first_model_estims <- extract_estims(boot_sample_model, beta_true, alpha)
    all_boot_betas[i] <- ifelse(nrow(first_model_estims) == 1, first_model_estims$beta_hat, NA)
    
    # Draws nboot_t number of samples from current bootstrap sample, fits model from nested sample,
    # and adds nested beta to separate list
    for (j in 1:nboot_t) {
      nested_boot_sample <- slice_sample(boot_sample, n = sample_size, replace = TRUE)
      nested_boot_sample_model <- fit_model(data = nested_boot_sample)
      nested_model_estims <- extract_estims(model = nested_boot_sample_model, 
                                            beta_true = beta_true, 
                                            alpha = alpha)
      
      all_nested_boot_betas[j] <- ifelse(nrow(nested_model_estims) == 1, nested_model_estims$beta_hat, NA)
      
    }
    
    # Gets SE of the bth estimate
    se_stars[i] = sd(all_nested_boot_betas,  na.rm = TRUE) 
  }
  
  # Returns bootstrapped betas and their SEs
  return (list(boot_betas = all_boot_betas, se_stars = se_stars))
}