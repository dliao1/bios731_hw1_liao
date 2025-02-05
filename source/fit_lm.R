library(tidyverse)
library(doParallel)
library(foreach)
library(tictoc)

# Fits linear regression model
fit_model <- function(data) {
  model <- lm(y ~ x, data = data)
  return (model)
}