library(tidyverse)
library(doParallel)
library(foreach)
library(dplyr)
library(here)

source(here("source", "utils.R"))
options(pillar.sigfig = 15)

# Ensure necessary directories exist
if (!dir.exists(here("results"))) {
  dir.create(here("results"))
}
if (!dir.exists(here("results", "sim_wald"))) {
  dir.create(here("results", "sim_wald"))
}
if (!dir.exists(here("results", "sim_boot"))) {
  dir.create(here("results", "sim_boot"))
}

# Cluster setup
cl <- makeCluster(4)  # Use 8 cores
registerDoParallel(cl)

# Simulation setup
mc_err <- 0.01
cover <- 0.95
alpha <- 1 - 0.95
n_sim <- round((cover * (1 - cover)) / mc_err^2)  

n <- c(10, 50, 500)
beta_true <- c(0, 0.5, 2)
err_type <- c(0, 1)  # 1 = normal, 0 = right-skewed

param_grid <- expand.grid(
  n = n,
  beta_true = beta_true,
  err_type = err_type
)

set.seed(2025)

# Seeds for each loop
seeds <- floor(runif(n_sim, 1, 10000))

# Parallelized Wald results
wald_results <- foreach(
  i = 1:nrow(param_grid),
  .packages = c("tibble", "dplyr", "tidyverse", "broom", "here", "doParallel", "foreach")
) %dopar% {
  params <- param_grid[i, ]
  
  # Run simulations for n_sim in parallel
  all_wald_estim <- foreach(
    j = 1:n_sim,
    .combine = rbind,  # i want this to stay a list of dfs
    .packages = c("tibble", "dplyr", "tidyverse", "broom", "here", "doParallel", "foreach")
  ) %do% {
    set.seed(seeds[j])
    simdata <- gen_data(n = params$n, beta_true = params$beta_true, err_type = params$err_type)
    extract_estims(simdata, params$beta_true, alpha) # could do separate return statement for readability
  }
  
  # Save individual parameter combination results
  filename <- paste0("scenario_", i, ".RDA")
  save(all_wald_estim, file = here("results", "sim_wald", filename))
  
  return(all_wald_estim)
}

save(wald_results, file = here("results", "sim_wald", "all_sim_wald_scenarios.RDA"))

# Parallelized Bootstrap results
boot_percent_results <- foreach(
  i = 1:nrow(param_grid),
  .packages = c("tibble", "dplyr", "tidyverse", "broom", "here", "doParallel", "foreach")
) %dopar% {
  params <- param_grid[i, ]
  
  # Run simulations for n_sim in parallel
  all_boot_percent_estim <- foreach(
    j = 1:n_sim,
    .combine = rbind,  
    .packages = c("tibble", "dplyr", "tidyverse", "broom", "here", "doParallel", "foreach")
  ) %do% {
    set.seed(seeds[j])
    simdata <- gen_data(n = params$n, beta_true = params$beta_true, err_type = params$err_type)
    extract_estim_boot_percent(simdata, params$beta_true, alpha) # separate return statement for readability?
  }
  
  # Save individual parameter combination results
  filename <- paste0("scenario_", i, ".RDA")
  save(all_boot_percent_estim, file = here("results", "sim_boot", filename))
  
  return(all_boot_percent_estim)
}

save(boot_percent_results, file = here("results", "sim_boot", "all_sim_boot_scenarios.RDA"))


stopCluster(cl)
