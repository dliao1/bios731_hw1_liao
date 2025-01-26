# R Setup

library(tidyverse)
library(doParallel)
library(foreach)
library(dplyr)
library(here)

source(here("source", "utils.R"))
options(pillar.sigfig = 15)

if (!dir.exists(here("results", "sim_wald"))) {
  dir.create(here("results", "sim_wald"))
}

# Cluster setup
cl <- makeCluster(4)  # Use 8 cores
registerDoParallel(cl)

# Simulation setup
mc_err <- 0.01
cover <- 0.95
alpha <- 1 - 0.95
n_sim <- round((cover * (1 - cover)) / mc_err^2)  # Round 


n = c(10, 50, 500)
beta_true = c(0, 0.5, 2)
err_type = c(0, 1) # 1 = normal, 0 = right skewed

param_grid = expand.grid(n = n,
                         n_sim = n_sim,
                         beta_true = beta_true,
                         err_type = err_type)


# Simulation setup

# Seeds for each loop
seeds = floor(runif(n_sim, 1, 10000))

# List of results where each element is a dataframe with the estimated betas, etc,
# for all 475 simulations (with n = 10/50/500) run for that set of parameters
# 
results <- foreach (i = 1:nrow(param_grid), .packages = c("tibble", "dplyr", "tidyverse", "broom")) %dopar% {
  params <- param_grid[i, ]
  all_sim_data <- vector("list", n_sim)
  for (j in 1:n_sim) {
    set.seed(seeds[j])
    simdata = gen_data(n = params$n,
                       beta_true = params$beta_true,
                       err_type = params$err_type)
    all_sim_data[[j]] <- extract_estim_wald(simdata, params$beta_true, alpha)
    
    
  }
  curr_params_data <- bind_rows(all_sim_data)
  filename = paste0("scenario_", i, ".RDA")
  save(curr_params_data,
       file = here::here("results", "sim_wald", filename))
  return (curr_params_data)
}
stopCluster(cl)



