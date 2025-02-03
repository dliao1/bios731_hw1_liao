library(tidyverse)
library(doParallel)
library(foreach)
library(dplyr)
library(here)
library(tictoc)

source(here("source", "gen_data.R"))
source(here("source", "fit_lm.R"))
source(here("source", "run_bootstraps.R"))
source(here("source", "get_estimates.R"))

options(pillar.sigfig = 15)

# Makes sure directories exist
if (!dir.exists(here("results"))) {
  dir.create(here("results"))
}
if (!dir.exists(here("results", "sim_wald"))) {
  dir.create(here("results", "sim_wald"))
}
if (!dir.exists(here("results", "sim_boot_percentile"))) {
  dir.create(here("results", "sim_boot_percentile"))
}
if (!dir.exists(here("results", "sim_boot_t"))) {
  dir.create(here("results", "sim_boot_t"))
}

if (!dir.exists(here("results", "sim_data"))) {
  dir.create(here("results", "sim_data"))
}


# Cluster setup
cl <- makeCluster(12)  # Use 12 cores
registerDoParallel(cl)

# Simulation setup
mc_err <- 0.01
cover <- 0.95
alpha <- 1 - 0.95
n_sim <- round((cover * (1 - cover)) / mc_err^2)

n <- c(10, 50, 500)
beta_true <- c(0, 0.5, 2)
err_type <- c(0, 1)  # 1 = normal, 0 = lognormal

param_grid <- expand.grid(
  n = n,
  beta_true = beta_true,
  err_type = err_type
)


# Seeds!
set.seed(3000)
seeds <- floor(runif(n_sim, 1, 10000))

# Iterates through 18 parameter combinations
for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  # Runs 475 simulations per parameter combo
  # sim_results structure: dataframe with 475 rows (simulations), and 4 columns (Wald/Bootstrap percentile/Bootstrap t/Simulation Data)
  # , each of which is a 1 row dataframe for the specified method (i.e. sim_result$wald[[10]] = beta hat and Wald CI for the 10th 
  # simulation run for the current parameter combination )

  # At the end of all 475 simulations for the current scenario, I convert every row
  # of sim_results to 1 row dataframes and bind them to dataframes for their corresponding method (Wald/percentile/etc.) 
  
  # End result: 
  # all_wald_estim = 1 dataframe, 475 rows
  
  sim_results <- foreach(
    j = 1:n_sim,
    .combine = rbind,
    .packages = c("tibble", "dplyr", "tidyverse", "broom", "here", "tictoc")
  ) %dopar% {
    set.seed(seeds[j])
    # Generates simulated data
    simdata <- gen_data(n = params$n, 
                        beta_true = params$beta_true, 
                        err_type = params$err_type
                        )
    
    # Computes Estimates
    # Note each: each result is one (1) single row
    
    model_fit <- fit_model(simdata)
    
    tic()
    wald_result <- extract_estims(model = model_fit, 
                                  beta_true = params$beta_true, 
                                  alpha = alpha)
    wald_time <- toc(quiet = TRUE)
    wald_result <- cbind(wald_result, scenario = i, sim = j, params, time = wald_time$toc - wald_time$tic)
    
    # Computes Bootstrap Percentile estimates
    nboot <- 50
    nboot_t <- 10
    
    tic()
    
    boot_data <- get_boot_data(original_data = simdata, 
                               beta_true = params$beta_true, 
                               sample_size = params$n,
                               nboot = nboot,
                               alpha = alpha)
    
    boot_percent_result <- extract_estim_boot_percent(all_boot_betas = boot_data,
                               beta_true = params$beta_true,
                               alpha = alpha)
    
    boot_percent_time <- toc(quiet = TRUE)
    
    boot_percent_result <- cbind(boot_percent_result, scenario = i, sim = j, params, time = boot_percent_time$toc - boot_percent_time$tic)
    
    # Boot t estimates
    tic()
    boot_t_data <- get_boot_t_data(original_data = simdata, 
                               beta_true = params$beta_true, 
                               sample_size = params$n,
                               alpha = alpha,
                               nboot = nboot,
                               nboot_t = nboot_t)
    
    boot_t_result <- extract_estim_boot_t(original_data = simdata,
                                          all_boot_betas = boot_t_data$boot_betas,
                                          se_stars = boot_t_data$se_stars,
                                          beta_true = params$beta_true,
                                          alpha = alpha)
    
    boot_t_time <- toc(quiet = TRUE)

    boot_t_result <- cbind(boot_t_result, scenario = i, sim = j, params, time = boot_t_time$toc - boot_t_time$tic)
    
    # Casts 2 rows into 2 lists and makes 2 columns, 1 for each list
    tibble(
      wald = list(wald_result),
      boot_percent = list(boot_percent_result),
      sim_data = list(simdata),
      boot_t = list(boot_t_result)
    )
  }
  
  # Turns each row in wald column into dataframe and binds all rows together for all 475 wald results for current
  #simulation
  all_wald_estim <- bind_rows(lapply(sim_results$wald, as.data.frame))
  all_boot_percent_estim <- bind_rows(lapply(sim_results$boot_percent, as.data.frame))
  all_sim_data <- bind_rows(lapply(sim_results$sim_data, as.data.frame))
  all_boot_t_estim <- bind_rows(lapply(sim_results$boot_t, as.data.frame))
  
  # Save **only** the single scenario’s results, not the full list!
  save(all_wald_estim, file = here("results", "sim_wald", paste0("scenario_", i, ".RDA")))
  save(all_boot_percent_estim, file = here("results", "sim_boot_percentile", paste0("scenario_", i, ".RDA")))
  save(all_sim_data, file = here("results", "sim_data", paste0("scenario_", i, ".RDA")))
  save(all_boot_t_estim, file = here("results", "sim_boot_t", paste0("scenario_", i, ".RDA")))
  
  
  # Print progress
  cat(sprintf("Saved scenario %d (n = %d, beta_true = %.2f, err_type = %d)\n",
              i, params$n, params$beta_true, params$err_type))
}

# Stop the cluster
stopCluster(cl)
