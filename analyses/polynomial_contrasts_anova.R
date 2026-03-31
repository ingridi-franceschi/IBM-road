############################################################
# STATISTICAL ANALYSIS — POLYNOMIAL RESPONSE TO TRAFFIC
#
# Description:
# This script evaluates linear, quadratic, and cubic trends in
# population growth across increasing traffic flow using
# orthogonal polynomial contrasts within an ANOVA framework.
#
# Note:
# Traffic flow is treated as an ordered factor with unequal spacing,
# and custom scores are used to reflect realistic traffic intensity.
############################################################

# ==========================================================
# 1. Load packages
# ==========================================================
library(dplyr)

# ==========================================================
# 2. Data aggregation
# ==========================================================
# Aggregate simulation outputs at the replicate level to avoid
# temporal pseudoreplication
analysis_data <- all_simulations_results |>
  group_by(replica, traffic_flow, cross_prob_class, category) |>
  summarise(
    mean_growth = mean(relative_growth, na.rm = TRUE),
    .groups = "drop"
  )

# ==========================================================
# 3. Define ordered traffic factor with polynomial contrasts
# ==========================================================
# Traffic levels are unevenly spaced; custom scores preserve
# the actual gradient used in simulations
traffic_levels <- c(0, 1, 5, 10)

analysis_data$traffic_flow <- factor(
  analysis_data$traffic_flow,
  levels = traffic_levels,
  ordered = TRUE
)

contrasts(analysis_data$traffic_flow) <- contr.poly(n = length(traffic_levels),
                                                    scores = traffic_levels)

# ==========================================================
# 4. Model 1 — Overall effect of traffic
# ==========================================================
# Tests whether population relative growth exhibits linear or non-linear
# responses to increasing traffic intensity
model_main <- aov(mean_growth ~ traffic_flow,
                  data = analysis_data)

summary.lm(model_main)

# ==========================================================
# 5. Model 2 — Interaction effects
# ==========================================================
# Evaluate whether traffic effects depend on:
# - crossing probability
# - life-history category
model_interaction <- aov(mean_growth ~ traffic_flow *
                           factor(cross_prob_class) *
                           factor(category),
                         data = analysis_data)

summary.lm(model_interaction)