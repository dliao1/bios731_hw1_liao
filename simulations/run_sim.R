# R Setup

library(tidyverse)
library(doParallel)
library(foreach)
library(dplyr)

source(here("source", "utils.R"))

# Cluster setup
cl <- makeCluster(4)  # Use 4 cores
registerDoParallel(cl)

# Simulation setup
mc_err <- 0.01
cover <- 0.95
alpha <- 1 - 0.95
n_sim <- (cover * (1 - cover))/mc_err^2


n = c(10, 50, 500)
beta_true = c(0, 0.5, 2)
err_type = c(0, 1) # 1 = normal, 0 = right skewed

param_grid = expand.grid(n = n,
                     n_sim = n_sim,
                     beta_true = beta_true,
                     err_type = err_type)

 # Delete when you run entire thing and replace with n_sim in for loop

# Simulation setup

# Seeds for each loop
seeds = floor(runif(n_sim, 1, 10000))

n_sim <- 2
param_grid <- param_grid[1, ]

# List of results where each element is all 750 simulations run for that set of parameters
results <- foreach (i = 1:nrow(param_grid), .packages = c("tibble", "dplyr")) %dopar% {
  set.seed(seeds[[i]])
  params <- param_grid[i, ]
  
  all_sim_data = as.list(rep(NA, n_sim))
  for (j in 1:n_sim) {
    simdata = gen_data(n = params$n,
                       beta_true = params$beta_true,
                       err_type = params$err_type)
    all_sim_data[[j]] <- extract_estim(simdata, params$beta_true)
    
    
  }
  curr_params_data <- bind_rows(all_sim_data)
  return (curr_params_data)
}
stopCluster(cl)



