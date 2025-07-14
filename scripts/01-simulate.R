library(tidyverse)
library(here)
library(truncnorm)

set.seed(123)

n <- 1000

simulated_patients <- tibble(
  age = sample(18:85, n, replace = TRUE),
  sex = sample(c("Male", "Female"), n, replace = TRUE),
  weight = round(rnorm(n, mean = 75, sd = 15), 1),       # kg
  height = round(rnorm(n, mean = 170, sd = 10), 1),      # cm
  creat_umol = round(runif(n, 40, 1000), 1)              # Creatinine in µmol/L
) %>%in
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
    crcl_lt_30 = crcl < 30
    )
  

grouped_gfr <- simulated_patients %>% 
  select(egfr, bsa, crcl, weight, delta) %>% 
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

ggplot(grouped_gfr, aes(x = bsa, y = delta, color = gfr_group)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 1.73, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  annotate("text", x = 1.4, y = 50, label = "CrCl > eGFR", hjust = 1, size = 4, alpha = 0.5) +
  annotate("text", x = 1.4, y = -50, label = "CrCl < eGFR", hjust = 1, size = 4, alpha = 0.5) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "green", alpha = 0.05) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "red", alpha = 0.05) +
    labs(title = "Effect of BSA on the difference between eGFR and CrCL",
       subtitle = "Delta = CrCL - eGFR",
       x = "BSA (m²)",
       y = "Delta",
       color = "eGFR Grouping (CKD-EPI)") +
  theme_minimal()

ggplot(simulated_patients, aes(x = egfr, fill = crcl_lt_30)) +
  geom_histogram(binwidth = 5, position = "stack", boundary = 0) +
  scale_fill_manual(values = c("grey70", "red")) +
  labs(
    title = "Distribution of eGFR and Risk of CrCl < 30",
    x = "eGFR (ml/min/1.73m²)",
    y = "Number of Patients",
    fill = "CrCl < 30"
  ) +
  theme_minimal()

thresholds <- seq(5, 90, by = 1)

threshold_summary <- thresholds %>%
  map_df(function(thresh) {
    data <- simulated_patients %>%
      mutate(flag = egfr < thresh)
    
    true_positives <- sum(data$flag & data$crcl_lt_30)
    false_negatives <- sum(!data$flag & data$crcl_lt_30)
    true_negatives <- sum(!data$flag & !data$crcl_lt_30)
    
    tibble(
      egfr_thresh = thresh,
      sensitivity = true_positives / (true_positives + false_negatives),
      npv = true_negatives / (true_negatives + false_negatives)
    )
  })

cutoffs <- threshold_summary %>% 
  filter(sensitivity >=0.95) %>% 
  summarise(
    egfr_95 = min(egfr_thresh[sensitivity >= 0.95]),
    egfr_99 = min(egfr_thresh[sensitivity >= 0.99]),
    egfr_100 = min(egfr_thresh[sensitivity == 1]),
  )

x_range <- range(threshold_summary$egfr_thresh, na.rm = TRUE)
y_range <- range(threshold_summary$sensitivity, na.rm = TRUE)

ggplot(threshold_summary, aes(x = egfr_thresh, y = sensitivity)) +
  geom_point(alpha = 0.6) +
  geom_line() +
  geom_segment(aes(x= cutoffs$egfr_95, xend = cutoffs$egfr_95, y= min(y_range), yend = 0.95),
                linetype = "dashed", color = "red") +
  geom_segment(aes(x= 0, xend = cutoffs$egfr_95, y= 0.95, yend = 0.95),
               linetype = "dashed", color = "red") +
  geom_segment(aes(x= cutoffs$egfr_99, xend = cutoffs$egfr_99, y= min(y_range), yend = 0.99),
               linetype = "dashed", color = "orange") +
  geom_segment(aes(x= 0, xend = cutoffs$egfr_99, y= 0.99, yend = 0.99),
               linetype = "dashed", color = "orange") +
  geom_segment(aes(x= cutoffs$egfr_100, xend = cutoffs$egfr_100, y= min(y_range), yend = 1),
               linetype = "dashed", color = "green") +
  geom_segment(aes(x= 0, xend = cutoffs$egfr_100, y= 1, yend = 1),
               linetype = "dashed", color = "green") +
  scale_x_continuous(breaks = seq(0, 90, by = 5), limits = c(0, NA)) +
  labs(
    title = "Sensitivity for CrCl > 30 for each eGFR threshold",
    x = "eGFR (ml/min/1.73m²)",
    y = "Sensitivity"
  ) +
  theme_classic()
