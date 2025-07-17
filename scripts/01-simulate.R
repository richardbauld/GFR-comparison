library(tidyverse)
library(here)
library(truncnorm)
library(ggpmisc)
library(janitor)

# data import
ww_export <- read_csv(here("data/age-bmi-export.csv"))

ww_export <- ww_export %>% 
  clean_names()

ww_summ <- ww_export %>%
  summarise(mean_age = mean(demographic_age_years),
            sd_age = sd(demographic_age_years),
            mean_height = mean(demographic_height, na.rm = TRUE),
            sd_height = sd(demographic_height, na.rm = TRUE),
            mean_weight = mean(demographic_weight, na.rm = TRUE),
            sd_weight = sd(demographic_weight, na.rm = TRUE))

set.seed(42)

n <- 1000

### Using data from HSE 2022 summary report ----
# N's
n_female <- 11527
n_male <- 10378
prop_male <- n_male / (n_male + n_female)
prop_female <- n_female / (n_male + n_female)

# Means
mean_h_male <- 176.2
mean_h_female <- 162.3
mean_w_male <- 85.8
mean_w_female <-72.8

# Standard errors
se_h_male <- 0.22
se_h_female <- 0.17
se_w_male <- 0.47
se_w_female <- 0.43

# SDs from SEs
sd_height_m <- se_h_male * sqrt(n_male)
sd_height_f <- se_h_female * sqrt(n_female)
sd_weight_m <- se_w_male * sqrt(n_male)
sd_weight_f <- se_w_female * sqrt(n_female)
  
  
simulated_patients <- tibble(
  age = rtruncnorm(n, a = 18, b = 90, mean = 58, sd = 16),
  sex = sample(c("Male", "Female"), size = n, replace = TRUE, prob = c(prop_male, prop_female)),
  weight = rtruncnorm(n, a = 30, b = 200,
                      mean = case_when(
                        sex == "Male" ~ mean_w_male,
                        sex == "Female" ~ mean_w_female
                      ),
                      sd = case_when(
                        sex == "Male" ~ sd_weight_m,
                        sex == "Female" ~ sd_weight_f
                      )),       # kg
  height = round(rtruncnorm(n,
                            a = case_when(
                              sex == "Male" ~ 155,
                              sex == "Female" ~ 145
                            ), 
                            b = case_when(
                              sex == "Male" ~ 205,
                              sex == "Female" ~ 190
                              ),
                       mean = case_when(
                         sex == "Male" ~ mean_h_male,
                         sex == "Female" ~ mean_h_female
                         ),
                       sd = case_when(
                         sex == "Male" ~ sd_height_m,
                         sex == "Female" ~ sd_height_f
                       )), 1),      # cm
  creat_umol = round(rtruncnorm(n, a = 40, b = 1200, mean = 120, sd = 400), 1)            # Creatinine in µmol/L
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
    ),
    delta = crcl - egfr,
    crcl_lt_30 = crcl < 30,
    corr_egfr = egfr * (bsa/1.73),
    mean_gfr = (crcl + egfr) / 2,
    flag = egfr > 30 & crcl < 30,
    corr_gfr_lt_30 = corr_egfr < 30
    )
  
mean_diff <- simulated_patients %>% 
  summarise(mean_diff = round(mean(delta), 2),
            sd_diff = sd(delta, na.rm = TRUE),
            lower_limit = mean_diff - 1.96 * sd_diff,
            upper_limit = mean_diff + 1.96 * sd_diff)

grouped_gfr <- simulated_patients %>% 
  select(egfr, bsa, crcl, weight, delta, flag) %>% 
  mutate(gfr_group = case_when(
    egfr < 30 ~ "<30",
    egfr >= 30 & egfr < 45 ~ "30–44",
    egfr >= 45 & egfr < 60 ~ "45–59",
    egfr >= 60 & egfr < 90 ~ "60–89",
    egfr >= 90 ~ "≥90"),
    gfr_group = factor(gfr_group, levels = c("<30", "30–44", "45–59", "60–89", "≥90")))

# Grab names of df's I want to keep
keep_names <- c("grouped_gfr", "mean_diff", "simulated_patients")

# Remove all objects except the data frame you want to keep
rm(list = setdiff(ls(), keep_names))


  









