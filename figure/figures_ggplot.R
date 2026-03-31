############################################################
# INDIVIDUAL-BASED MODEL (IBM)
# Behavioural and life-history traits and road traffic drive
# population persistence: a mechanistic modelling framework
#
# Date: 2026-03-31
# Version: 1.0
#
# Description:
# This script generates Figures 1–3 summarizing how population
# growth and collision risk respond to traffic flow, crossing
# probability, and life-history strategies.
#
# Packages: dplyr, ggplot2, ggdist, ggh4x, scales
############################################################

# ==========================================================
# 1. Load required packages
# ==========================================================
library(dplyr)
library(ggplot2)
library(ggdist)
library(ggh4x)
library(scales)

# ==========================================================
# 2. Helper variables (life-history classification)
# ==========================================================
# Define life-history groupings used across all figures
data <- all_simulations_results |>
  mutate(
    mobility_class = ifelse(category %in% c(1, 4),
                            "Low daily movement",
                            "Large daily movement"),
    mature_class = ifelse(category %in% c(1, 2),
                          "Early mature",
                          "Late mature")
  )

# ==========================================================
# 3. FIGURE 1 — Growth response to traffic and crossing probability
# ==========================================================
# Summarize post-transient dynamics (year >= 6)
fig1_data <- data |>
  filter(year >= 6) |>
  group_by(mature_class, mobility_class,
           traffic_flow, cross_prob_class, replica) |>
  summarise(
    collision_rate = mean(collision_rate, na.rm = TRUE),
    mean_growth = mean(relative_growth, na.rm = TRUE),
    mean_birth = mean(birth_total, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    cross_prob_class = factor(cross_prob_class,
                              levels = c("low", "medium", "high"))
  )

fig1 <- ggplot(fig1_data,
               aes(x = factor(traffic_flow),
                   y = mean_growth,
                   color = cross_prob_class,
                   group = cross_prob_class)) +
  
  # Smoothed trends highlight non-linear responses to traffic
  geom_smooth(size = 1, se = FALSE) +
  
  # Mean values per traffic level
  stat_summary(fun = mean, geom = "point", size = 3) +
  
  # Reference line for population stability
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  facet_grid(mature_class ~ mobility_class, scales = "free_y") +
  
  scale_color_manual(values = c(
    "low" = "#7570b3",
    "medium" = "#018571",
    "high" = "#f1a340"
  )) +
  
  labs(x = "Traffic flow (vehicles/min)",
       y = "Mean relative growth",
       color = "Crossing probability") +
  
  theme_minimal(base_size = 25) +
  theme(strip.text = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.spacing = unit(1, "cm"),
        axis.ticks = element_line(color = "black"))

# Apply facet-specific y-axis scaling
fig1 <- fig1 + ggh4x::facetted_pos_scales(
  y = list(
    mature_class == "Early mature" ~
      scale_y_continuous(labels = label_number(drop0trailing = TRUE),
                         limits = c(-0.0005, 0.03)),
    mature_class == "Late mature" ~
      scale_y_continuous(labels = label_number(drop0trailing = TRUE),
                         limits = c(-0.01, 0.015))
  )
)

ggsave("fig01.png", fig1,
       width = 50, height = 20, units = "cm",
       dpi = 300, bg = "white")

# ==========================================================
# 4. FIGURE 2 — Temporal dynamics of population growth
# ==========================================================
# Evaluate how growth trajectories evolve through time
fig2_data <- data |>
  group_by(mobility_class, mature_class,
           traffic_flow, year) |>
  summarise(
    mean_growth = mean(relative_growth, na.rm = TRUE),
    .groups = "drop"
  )

fig2 <- ggplot(fig2_data,
               aes(x = year,
                   y = mean_growth,
                   color = factor(traffic_flow),
                   group = factor(traffic_flow))) +
  
  geom_line(size = 1) +
  
  # Reference lines: equilibrium and transient cutoff
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 6, linetype = "dashed") +
  
  facet_grid(mature_class ~ mobility_class, scales = "free_y") +
  
  scale_color_manual(values = c(
    "0" = "black",
    "1" = "#7570b3",
    "5" = "#018571",
    "10" = "#f1a340"
  )) +
  
  labs(x = "Years",
       y = "Mean relative growth",
       color = "Traffic flow (vehicles/min)") +
  
  scale_x_continuous(limits = c(1, 20),
                     breaks = c(1, 5, 10, 15, 20),
                     expand = c(0, 0)) +
  
  theme_minimal(base_size = 25) +
  theme(strip.text = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA),
        panel.spacing = unit(1, "cm"),
        axis.ticks = element_line(color = "black"))

# Facet-specific scaling improves interpretability across strategies
fig2 <- fig2 + ggh4x::facetted_pos_scales(
  y = list(
    mature_class == "Early mature" ~
      scale_y_continuous(limits = c(-0.2, 0.6),
                         breaks = seq(-0.2, 0.6, 0.2)),
    mature_class == "Late mature" ~
      scale_y_continuous(limits = c(-0.1, 0.1),
                         breaks = seq(-0.1, 0.1, 0.05))
  )
)

ggsave("fig02.png", fig2,
       width = 50, height = 20, units = "cm",
       dpi = 300, bg = "white")

# ==========================================================
# 5. FIGURE 3 — Collision risk distribution
# ==========================================================
# Focus on scenarios with traffic flow
fig3_data <- data |>
  filter(year >= 6, traffic_flow != 0) |>
  group_by(mature_class, mobility_class,
           traffic_flow, cross_prob_class, replica) |>
  summarise(
    collision_rate = mean(collision_rate, na.rm = TRUE),
    mean_growth = mean(relative_growth, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    cross_prob_class = factor(cross_prob_class,
                              levels = c("low", "medium", "high"))
  )

fig3 <- ggplot(fig3_data,
               aes(x = cross_prob_class,
                   y = collision_rate,
                   fill = factor(traffic_flow),
                   color = factor(traffic_flow))) +
  
  # Distribution + uncertainty (posterior-like visualization)
  stat_halfeye(
    justification = 0,
    point_interval = mean_qi,
    position = position_dodge(width = 0.6),
    alpha = 0.5
  ) +
  
  # Mean and interval estimates
  stat_pointinterval(
    point_interval = mean_qi,
    position = position_dodge(width = 0.6),
    .width = 0.95,
    size = 3
  ) +
  
  facet_grid(mature_class ~ mobility_class) +
  
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  scale_fill_manual(values = c(
    "1" = "#7570b3",
    "5" = "#018571",
    "10" = "#f1a340"
  )) +
  
  scale_color_manual(values = c(
    "1" = "#7570b3",
    "5" = "#018571",
    "10" = "#f1a340"
  )) +
  
  guides(color = "none") +
  
  labs(x = "Crossing probability",
       y = "Mean daily collision rate",
       fill = "Traffic flow (vehicles/min)") +
  
  theme_minimal(base_size = 25) +
  theme(strip.text = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        panel.border = element_rect(color = "black", fill = NA),
        axis.ticks = element_line(color = "black"))

ggsave("fig03.png", fig3,
       width = 50, height = 20, units = "cm",
       dpi = 300, bg = "white")