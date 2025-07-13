library(tidyverse)
library(here)

set.seed(123)

n <- 500

simulated_patients <- tibble(
  age = sample(18:85, n, replace = TRUE),
  sex = sample(c("Male", "Female"), n, replace = TRUE),
  weight = round(rnorm(n, mean = 75, sd = 15), 1),       # kg
  height = round(rnorm(n, mean = 170, sd = 10), 1),      # cm
  creat_umol = round(runif(n, 40, 1000), 1)              # Creatinine in µmol/L
) %>%
  mutate(
    creat_mgdl = creat_umol / 88.4,
    kappa = if_else(sex == "Female", 0.7, 0.9),
    alpha = if_else(sex == "Female", -0.241, -0.302),
    sex_coeff = if_else(sex == "Female", 1.012, 1),
    cr_over_kappa = creat_mgdl / kappa,
    egfr = 142 * 
      (pmin(cr_over_kappa, 1) ^ alpha) *
      (pmax(cr_over_kappa, 1) ^ (-1.200)) *
      (0.9938 ^ age) *
      sex_coeff,
    crcl = if_else(
      sex == "Male",
      ((140 - age) * weight) / (72 * creat_mgdl),
      0.85 * ((140 - age) * weight) / (72 * creat_mgdl)
    )
  )
