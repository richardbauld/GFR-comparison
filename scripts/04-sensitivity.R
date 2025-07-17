# Run this if not already run
source(here("scripts/01-simulate.R"))


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
