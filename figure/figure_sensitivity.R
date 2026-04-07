############################################################
# SOBOL SENSITIVITY ANALYSIS — FIGURES
#
# Description:
# This script generates Sobol sensitivity plots including:
# (1) First- (Si) and total-order (Ti) indices
# (2) Second-order interaction indices (Sij)
#
# These results quantify the relative importance of life-
# history traits and traffic-related parameters on model output.
############################################################

# ==========================================================
# 1. Load packages
# ==========================================================
library(dplyr)
library(ggplot2)
library(tidyr)

# ==========================================================
# 2. Parameter labels
# ==========================================================
param_labels <- c(
  "lmin"           = "Min. step length",
  "lmax"           = "Max. step length",
  "mature"         = "Age at maturity",
  "num_offspring"  = "No. offspring",
  "K"              = "Carrying capacity",
  "min_prob_cross" = "Min. crossing prob.",
  "max_prob_cross" = "Max. crossing prob.",
  "traffic_flow"   = "Traffic flow"
)

# Define plotting order (importance-driven layout)
param_order <- c(
  "Traffic flow",
  "Max. crossing prob.",
  "Min. crossing prob.",
  "Carrying capacity",
  "No. offspring",
  "Age at maturity",
  "Max. step length",
  "Min. step length"
)

# ==========================================================
# 3. FIGURE 4 — First- and total-order Sobol indices
# ==========================================================
# Filter and relabel sensitivity indices
fig4_data <- indice_sestiv$results |>
  filter(sensitivity %in% c("Si", "Ti")) |>
  mutate(
    parameters = recode(parameters, !!!param_labels),
    parameters = factor(parameters, levels = param_order)
  )

fig4 <- ggplot(fig4_data,
               aes(x = parameters, y = original)) +
  
  # Bar plot of Sobol indices
  geom_col(fill = "steelblue", width = 0.6) +
  
  # Confidence intervals (uncertainty from Sobol estimation)
  geom_errorbar(aes(ymin = low.ci, ymax = high.ci),
                width = 0.3,
                linewidth = 1,
                color = "black") +
  
  # Reference line (no effect)
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  
  # Separate panels for Si and Ti
  facet_wrap(~ sensitivity,
             labeller = as_labeller(c(
               Si = "First-order (Sᵢ)",
               Ti = "Total-order (Tᵢ)"
             ))) +
  
  coord_flip() +
  
  labs(x = NULL,
       y = "Sobol' index") +
  
  theme_bw(base_size = 40) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

ggsave("fig04_sobol_main.png", fig4,
       width = 30, height = 20, units = "cm",
       dpi = 300, bg = "white")

# ==========================================================
# 4. FIGURE 5 — Second-order interaction effects (Sij)
# ==========================================================
# Extract pairwise parameter interactions
fig5_data <- indice_sestiv$results |>
  filter(sensitivity == "Sij") |>
  separate(parameters, into = c("param1", "param2"), sep = "\\.") |>
  mutate(
    param1 = recode(param1, !!!param_labels),
    param2 = recode(param2, !!!param_labels),
    
    # Create readable interaction labels
    pair = paste(param1, "×", param2),
    
    facet_label = "Second-order (Sᵢⱼ)"
  )

fig5 <- ggplot(fig5_data,
               aes(x = reorder(pair, original),
                   y = original)) +
  
  # Interaction strength
  geom_col(fill = "steelblue", width = 0.6) +
  
  # Uncertainty intervals
  geom_errorbar(aes(ymin = low.ci, ymax = high.ci),
                width = 0.3,
                linewidth = 1,
                color = "black") +
  
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "grey50") +
  
  facet_wrap(~ facet_label) +
  
  coord_flip() +
  
  labs(x = NULL,
       y = "Sobol' index") +
  
  theme_bw(base_size = 40) +
  theme(
    strip.text = element_text(face = "bold"),
    panel.border = element_rect(color = "black", fill = NA),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

ggsave("fig05_sobol_interactions.png", fig5,
       width = 30, height = 25, units = "cm",
       dpi = 300, bg = "white")