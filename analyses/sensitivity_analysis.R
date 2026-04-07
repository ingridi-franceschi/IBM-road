############################################################
# SOBOL SENSITIVITY ANALYSIS — MODEL RUN
#
# Description:
# This script performs a global sensitivity analysis using
# Sobol variance-based indices. It:
# (1) Generates Sobol sampling matrices
# (2) Scales parameters to biologically realistic ranges
# (3) Runs the IBM model for each parameter set
# (4) Computes first-, total-, and second-order indices
#
# Output:
# - CSV file with simulation results
# - Sobol sensitivity indices object
############################################################

# ==========================================================
# 1. Load packages
# ==========================================================
library(sensobol)
library(data.table)
library(dplyr)
library(tidyr)
library(filelock)

# ==========================================================
# 2. Define variable parameters
# ==========================================================
params_variaveis <- c(
  "lmin", "lmax", "mature", "K",
  "num_offspring", "traffic_flow",
  "min_prob_cross", "max_prob_cross"
)

# ==========================================================
# 3. Generate Sobol sampling matrix
# ==========================================================
N <- 100  # samples per base matrix

# Total rows:
# L = N * (2 + P + P(P-1)/2)  (second-order design)

mat <- sensobol::sobol_matrices(
  matrices = c("A", "B", "AB"),
  N = N,
  params = params_variaveis,
  order = "second"
)

# ==========================================================
# 4. Scale parameters to biological ranges
# ==========================================================

# Movement
mat[, "lmin"] <- runif(mat[, "lmin"], 6.94, 69.44)
mat[, "lmax"] <- pmax(mat[, "lmin"] + 0.05, runif(nrow(mat), 34.72, 104.17))

# Life-history
mat[, "mature"] <- round(runif(mat[, "mature"], 1,3))
mat[, "K"] <- round(runif(mat[, "K"], 50,100))
mat[, "num_offspring"] <- round(runif(mat[, "num_offspring"], 1, 8))

# Traffic
mat[, "traffic_flow"] <- round(runif(mat[, "traffic_flow"], 1,10))

# Behaviour
mat[, "min_prob_cross"] <- runif(mat[, "min_prob_cross"], 0, 0.66)
mat[, "max_prob_cross"] <- pmax(mat[, "min_prob_cross"] + 0.05, runif(nrow(mat),0.33, 1))

# ==========================================================
# 5. Fixed parameters
# ==========================================================
params_fixos <- matrix(
  c(
    144, 1, 0, pi/2,
    1, 1, 84, 100
  ),
  nrow = nrow(mat),
  ncol = 8,
  byrow = TRUE
)

colnames(params_fixos) <- c(
  "steps", "mu", "angle_mean", "angle_sd",
  "reprod_period", "reprod_freq",
  "speed_mean", "speed_sd"
)

# ==========================================================
# 6. Combine parameters
# ==========================================================
params <- cbind(mat, params_fixos)

all_pars <- as.data.frame(params) |>
  as_tibble() |>
  relocate("steps", .before = "lmin") |>
  relocate("mu", .after = "steps") |>
  relocate("angle_mean", .after = "lmax") |>
  relocate("angle_sd", .after = "angle_mean") |>
  relocate("reprod_period", .after = "mature") |>
  relocate("reprod_freq", .after = "reprod_period")

# ==========================================================
# 7. Run IBM model
# ==========================================================
source("./IBM_models_sensibility.R")

if (!dir.exists("output")) dir.create("output")
if (!dir.exists("logs")) dir.create("logs")

args <- commandArgs(trailingOnly = TRUE)
idx <- as.integer(args[1])

# Files
lock_csv <- "./output/all_results_sens_final1.lock"
csv_file <- "./output/all_results_sens_final1.csv"

# Initialize output file
lock <- lock(lock_csv)
if (!file.exists(csv_file)) {
  empty_dt <- data.table(
    sim = integer(),
    year = integer(),
    initial_population = numeric(),
    end_population = numeric(),
    crossings = numeric(),
    crossing_rate = numeric(),
    collision_rate = numeric(),
    death_collision = numeric(),
    death_random = numeric(),
    birth_total = numeric(),
    net_growth = numeric(),
    relative_growth = numeric()
  )
  fwrite(empty_dt, csv_file)
}
unlock(lock)

# Logging
log_file <- "./output/simulation.txt"
lock_log <- "./output/simulation.lock"

lock_l <- lock(lock_log)
write(paste0("Simulation ", idx, " started at ", Sys.time()),
      file = log_file, append = TRUE)
unlock(lock_l)

# Run model
row <- all_pars[idx, ]

result <- IBM(
  b = 0.5,
  d = 0.2,
  K = row$K,
  N0 = 50,
  days = 365,
  years = 5,
  nrow = 2000,
  ncol = 2000,
  res = 50,
  road_width = 5,
  pars = row,
  sim = idx,
  lock_csv, csv_file, log_file, lock_log
)

# Log end
lock_l <- lock(lock_log)
write(paste0("Simulation ", idx, " completed at ", Sys.time()),
      file = log_file, append = TRUE)
unlock(lock_l)

# ==========================================================
# 8. Aggregate model output
# ==========================================================
data_sensit <- read.csv("all_results_sens_final1.csv")

resu <- data_sensit |>
  group_by(sim) |>
  summarise(mean_relative_growth = mean(relative_growth))

y <- as.vector(resu$mean_relative_growth)

# ==========================================================
# 9. Sobol sensitivity indices
# ==========================================================
indice_sestiv <- sensobol::sobol_indices(
  matrices = c("A", "B", "AB"),
  Y = y,
  N = N,
  params = params_variaveis,
  boot = TRUE,
  R = 1e5,
  type = "norm",
  first = "saltelli",
  total = "sobol",
  order = "second"
)