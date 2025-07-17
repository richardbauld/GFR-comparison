library(tidyverse)
library(here)
library(truncnorm)
library(ggpmisc)

set.seed(42)

n <- 1000

simulated_patients <- tibble(
  age = rtruncnorm(n, a = 18, b = 90, mean = 60, sd = 10),
  sex = sample(c("Male", "Female"), size = n, replace = TRUE, prob = c(0.45, 0.55)),
  weight = rtruncnorm(n, a = 30, b = 150, mean = 70, sd = 12),       # kg
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


ggplot(simulated_patients, aes(x = egfr, y = crcl, color = sex)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Comparison of eGFR (CKD-EPI 2021) and CrCl (Cockcroft-Gault)",
       x = "eGFR (mL/min/1.73m²)",
       y = "CrCl (mL/min)") +
  theme_minimal()

ggplot(simulated_patients, aes(x = corr_egfr, y = crcl, color = sex)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Comparison of eGFR corrected with BSA (CKD-EPI 2021) and CrCl (Cockcroft-Gault)",
       x = "Corrected eGFR (mL/min)",
       y = "CrCl (mL/min)") +
  theme_minimal()

ggplot(grouped_gfr, aes(x = bsa, y = crcl, color = gfr_group)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 30, linetype = "dashed", color = "black") +
  annotate("text", x = 2.25, y = 35, label = "CrCl = 30", hjust = 1, size = 4, alpha = 0.5) +
  labs(title = "Comparison of CrCl (Cockcroft-Gault) with BSA (DuBois formula)",
       subtitle = "By eGFR group (CKD-EPI)",
       x = "BSA (m²)",
       y = "CrCl (mL/min)",
       color = "GFR Grouping (CKD-EPI)") +
  theme_minimal()

ggplot(grouped_gfr, aes(x = bsa, y = delta, color = gfr_group, shape = flag)) +
  geom_point(alpha = 0.6) +
  geom_vline(xintercept = 1.73, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  annotate("text", x = 1.5, y = 50, label = "CrCl > eGFR", hjust = 1, size = 4, alpha = 0.5) +
  annotate("text", x = 1.5, y = -50, label = "CrCl < eGFR", hjust = 1, size = 4, alpha = 0.5) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "green", alpha = 0.05) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "red", alpha = 0.05) +
    labs(title = "Effect of BSA on the difference between eGFR and CrCL",
       subtitle = "Delta = CrCL - eGFR",
       x = "BSA (m²)",
       y = "Delta",
       color = "eGFR Grouping (CKD-EPI)",
       shape = "eGFR > 30 and CrCl < 30") +
  theme_minimal()

#ungrouped
ggplot(grouped_gfr, aes(x = bsa, y = delta)) +
  geom_point(alpha = 0.6) +
  geom_smooth() +
  geom_vline(xintercept = 1.73, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  annotate("text", x = 1.5, y = 50, label = "CrCl > eGFR", hjust = 1, size = 4, alpha = 0.5) +
  annotate("text", x = 1.5, y = -50, label = "CrCl < eGFR", hjust = 1, size = 4, alpha = 0.5) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf,
           fill = "green", alpha = 0.05) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
           fill = "red", alpha = 0.05) +
  stat_poly_eq(
    aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
    formula = y ~ x,
    method = "lm",
    parse = TRUE,
    npcx = 0.95,
    npcy = 0.5
  ) +
  labs(title = "Effect of BSA on the difference between eGFR and CrCL",
       subtitle = "Delta = CrCL - eGFR",
       x = "BSA (m²)",
       y = "Delta",
       color = "eGFR Grouping (CKD-EPI)",
       shape = "eGFR > 30 and CrCl < 30") +
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

# Need to create gfr thresholds to later test
thresholds <- seq(5, 90, by = 1)

# Function to run the gfr through each threshold and flag
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

# Creates cut-offs for each threshold to plot
cutoffs <- threshold_summary %>% 
  filter(sensitivity >=0.95) %>% 
  summarise(
    egfr_95 = min(egfr_thresh[sensitivity >= 0.95]),
    egfr_99 = min(egfr_thresh[sensitivity >= 0.99]),
    egfr_100 = min(egfr_thresh[sensitivity == 1]),
  )

# Generate data ranges for segment plots (so lines don't go further than the curve)
x_range <- range(threshold_summary$egfr_thresh, na.rm = TRUE)
y_range <- range(threshold_summary$sensitivity, na.rm = TRUE)

# Sensitivity plot
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
  theme_minimal()



# Bland-altman plot
ggplot(simulated_patients, aes(x = mean_gfr, y = delta)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = mean_diff$mean_diff, linetype = "dashed") +
  geom_hline(yintercept = mean_diff$lower_limit, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = mean_diff$upper_limit, linetype = "dashed", color = "orange") +
  annotate("text", x = 150, y = mean_diff$upper_limit, label = "+ 1.96 SD", hjust = 1, size = 4, alpha = 0.5, vjust = -1) +
  annotate("text", x = 150, y = mean_diff$lower_limit, label = "- 1.96 SD", hjust = 1, size = 4, alpha = 0.5, vjust = 1.5) +
  annotate("text", x = 150, y = mean_diff$mean_diff, label = "Mean", hjust = 1, size = 4, alpha = 0.5, vjust = -1) +
  annotate("text", x = 150, y = mean_diff$mean_diff, label = mean_diff$mean_diff, hjust = 1, size = 4, alpha = 0.5, vjust = 1.5) +
  labs(title = "Bland-Altman Plot: CrCl vs eGFR",
       x = "Mean of CrCl and eGFR",
       y = "Difference (CrCl − eGFR)") +
  theme_minimal()

ggplot(simulated_patients, aes(x = age)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
  labs(title = "Distribution of CrCl − eGFR",
       x = "Delta",
       y = "Number of Patients") +
  theme_minimal()

ggplot(simulated_patients, aes(x = height)) +
  geom_histogram(binwidth = 5, fill = "steelblue", color = "white") +
  labs(title = "Distribution of CrCl − eGFR",
       x = "Delta",
       y = "Number of Patients") +
  theme_minimal()

ggplot(simulated_patients, aes(x = age, y = delta)) + 
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm") +
  stat_poly_eq(
    aes(label = paste(after_stat(rr.label), after_stat(p.value.label), sep = "~~~")),
    formula = y ~ x,
    method = "lm",
    parse = TRUE,
    npcx = 0.95,
    npcy = 0.5
  ) +
  labs(title = "Delta depending on age ",
       x = "Age (Years)",
       y = "Delta") +
  theme_minimal()

# Sensitivity using eGFR corrected to bsa
# Function to run the gfr through each threshold and flag
threshold_summary_corr <- thresholds %>%
  map_df(function(thresh) {
    data <- simulated_patients %>%
      mutate(flag = corr_egfr < thresh)
    
    true_positives <- sum(data$flag & data$crcl_lt_30)
    false_negatives <- sum(!data$flag & data$crcl_lt_30)
    true_negatives <- sum(!data$flag & !data$crcl_lt_30)
    
    tibble(
      egfr_thresh = thresh,
      sensitivity = true_positives / (true_positives + false_negatives),
      npv = true_negatives / (true_negatives + false_negatives)
    )
  })

# Creates cut-offs for each threshold to plot
cutoffs_corr <- threshold_summary_corr %>% 
  filter(sensitivity >=0.95) %>% 
  summarise(
    egfr_95 = min(egfr_thresh[sensitivity >= 0.95]),
    egfr_99 = min(egfr_thresh[sensitivity >= 0.99]),
    egfr_100 = min(egfr_thresh[sensitivity == 1]),
  )

# Generate data ranges for segment plots (so lines don't go further than the curve)
x_range_corr <- range(threshold_summary_corr$egfr_thresh, na.rm = TRUE)
y_range_corr <- range(threshold_summary_corr$sensitivity, na.rm = TRUE)

# Sensitivity plot
ggplot(threshold_summary_corr, aes(x = egfr_thresh, y = sensitivity)) +
  geom_point(alpha = 0.6) +
  geom_line() +
  geom_segment(aes(x= cutoffs_corr$egfr_95, xend = cutoffs_corr$egfr_95, y= min(y_range_corr), yend = 0.95),
               linetype = "dashed", color = "red") +
  geom_segment(aes(x= 0, xend = cutoffs_corr$egfr_95, y= 0.95, yend = 0.95),
               linetype = "dashed", color = "red") +
  geom_segment(aes(x= cutoffs_corr$egfr_99, xend = cutoffs_corr$egfr_99, y= min(y_range_corr), yend = 0.99),
               linetype = "dashed", color = "orange") +
  geom_segment(aes(x= 0, xend = cutoffs_corr$egfr_99, y= 0.99, yend = 0.99),
               linetype = "dashed", color = "orange") +
  geom_segment(aes(x= cutoffs_corr$egfr_100, xend = cutoffs_corr$egfr_100, y= min(y_range_corr), yend = 1),
               linetype = "dashed", color = "green") +
  geom_segment(aes(x= 0, xend = cutoffs_corr$egfr_100, y= 1, yend = 1),
               linetype = "dashed", color = "green") +
  scale_x_continuous(breaks = seq(0, 90, by = 5), limits = c(0, NA)) +
  labs(
    title = "Sensitivity for CrCl > 30 for each BSA corrected eGFR threshold",
    x = "BSA corrected eGFR (ml/min)",
    y = "Sensitivity"
  ) +
  theme_minimal()
