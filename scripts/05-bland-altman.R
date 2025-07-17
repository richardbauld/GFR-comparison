# Run this if not already run
source(here("scripts/01-simulate.R"))

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