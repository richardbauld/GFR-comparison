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
    bsa = 0.007184 * (height^0.725) * (weight^0.425),
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

grouped_gfr <- simulated_patients %>% 
  select(egfr, bsa, crcl) %>% 
  mutate(gfr_group = case_when(
    egfr < 30 ~ "<30",
    egfr >= 30 & egfr < 45 ~ "30–44",
    egfr >= 45 & egfr < 60 ~ "45–59",
    egfr >= 60 & egfr < 90 ~ "60–89",
    egfr >= 90 ~ "≥90"),
    gfr_group = factor(gfr_group, levels = c("<30", "30–44", "45–59", "60–89", "≥90")))


ggplot(simulated_patients, aes(x = egfr, y = crcl, color = sex)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Comparison of eGFR (CKD-EPI 2021) and CrCl (Cockcroft-Gault)",
       x = "eGFR (mL/min/1.73m²)",
       y = "CrCl (mL/min)") +
  theme_minimal()

ggplot(grouped_gfr, aes(x = bsa, y = crcl, color = gfr_group)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 30, linetype = "dashed", color = "black") +
  labs(title = "Comparison of CrCl (Cockcroft-Gault) with BSA (DuBois formula)",
       subtitle = "By eGFR group (CKD-EPI)",
       x = "BSA (m²)",
       y = "CrCl (mL/min)",
       color = "GFR Grouping (CKD-EPI)") +
  theme_minimal()
