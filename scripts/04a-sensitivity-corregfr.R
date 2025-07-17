# Run this if not already run
source(here("scripts/01-simulate.R"))

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

