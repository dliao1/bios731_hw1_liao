library(tidyverse)


# test
test_func <- function(x) {
  result <- x^2
  return(result) 
}

gen_data <- function(n, beta_true, err_type) {
  beta0 <- 0
  beta_treat <- beta_true
  x <- rbinom(n, 1, prob = 0.5)
  epsilon <- ifelse (err_type == 1, rnorm(n, sd = sqrt(2)), rlnorm(10, meanlog = 0, sdlog = log(2)))
  y = beta0 + beta_treat * x + epsilon
  
  tibble(
    x = x,
    y = y
  )
}

extract_estim <- function(simdata, beta_true) {
  x <- simdata$x
  y <- simdata$y
  model <- lm(y ~ x)
  model_summary <- summary(model)
  
  beta_estim <- coef(model)["x"]  
  std_error <- model_summary$coefficients["x", "Std. Error"]
  conf_int <- confint(model)["x", ]
  
  coverage <- (ifelse(conf_int[1] <= beta_true & beta_true <= conf_int[2], 1, 0))
  
  tibble(
    beta_estim = beta_estim,
    se_beta = std_error,
    coverage = coverage
  )
}