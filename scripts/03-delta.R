# Run this if not already run
source(here("scripts/01-simulate.R"))

### Delta plot for BSA ----
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

### Ungrouped delta ----
# Includes regression line, etc
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

### Delta by age ----
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
