# R Setup

library(tidyverse)
library(doParallel)
library(foreach)
library(dplyr)
library(here)

source(here("source", "utils.R"))
options(pillar.sigfig = 15)

if (!dir.exists(here("results"))) {
  dir.create(here("results"))
}

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