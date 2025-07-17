# Run this if not already run
source(here("scripts/01-simulate.R"))

### CrCl vs GFR slope with sex groups ----
ggplot(simulated_patients, aes(x = egfr, y = crcl, color = sex)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Comparison of eGFR (CKD-EPI 2021) and CrCl (Cockcroft-Gault)",
       x = "eGFR (mL/min/1.73m²)",
       y = "CrCl (mL/min)") +
  theme_minimal()

### BSA corrected CKD-EPI ----
# Same as above but using BSA correction to remove indexing
ggplot(simulated_patients, aes(x = corr_egfr, y = crcl, color = sex)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(title = "Comparison of eGFR corrected with BSA (CKD-EPI 2021) and CrCl (Cockcroft-Gault)",
       x = "Corrected eGFR (mL/min)",
       y = "CrCl (mL/min)") +
  theme_minimal()

### CrCl vs BSA ----
# Using GFR groups
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

