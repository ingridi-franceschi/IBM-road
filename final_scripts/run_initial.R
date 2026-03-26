#------------------------------------------------------------------------------------------------------#
#------------------------- INDIVIDUAL-BASED MODEL (IBM) -----------------------------------------------#
#------------------------------------------------------------------------------#
# Ecological traits and road traffic drive population persistence: a mechanistic modelling framework
#
# Date: 2026
# Version: 1.0
#------------------------------------------------------------------------------------------------------#
# DESCRIPTION:
#   This script sets up and executes the IBM simulation across all parameter
#   combinations defined in the trait-by-traffic parameter table (pars).
#   It configures the landscape, population traits, road-crossing probability,
#   and traffic scenarios, then launches parallel simulations using the future.apply framework.
#
#   The parameter table (pars) contains 120 combinations derived from:
#     - 4 population categories (varying movement distance and life history)
#     - 3 road-crossing probability ranges (low, medium, high)
#     - 4 traffic flow levels (0, 1, 5, 10 vehicles/min)
#------------------------------------------------------------------------------------------------------#
# POPULATION CATEGORIES:
#   Four theoretical population types are defined based on combinations of
#   movement range and life history strategy:
#
#   Category 1: Short movement (1–5 km/day),  early maturity, high fecundity
#   Category 2: Long movement  (5–10 km/day), early maturity, high fecundity
#   Category 3: Long movement  (5–10 km/day), late maturity, low fecundity
#   Category 4: Short movement (1–5 km/day),  late maturity, low fecundity
#------------------------------------------------------------------------------------------------------#
# DEPENDENCIES:
#   tictoc       : Execution timing
#   data.table   : Fast CSV writing and reading
#   filelock     : File locking for safe parallel output writing
#   dplyr        : Parameter table manipulation
#   tidyr        : Parameter combination expansion (crossing())
#   future       : Parallel processing backend
#   future.apply : Parallel lapply interface
#
#   Install with:
#   install.packages(c("tictoc", "data.table", "filelock",
#                      "dplyr", "tidyr", "future", "future.apply"))
#------------------------------------------------------------------------------------------------------#
# ---- Load required packages ----
library(tictoc)
library(data.table)
library(filelock)

# ---- Load IBM model functions ----
source("./IBM_models_v0.R")

# ---- Create output directories if they do not exist ----
if (!dir.exists("output")) dir.create("output")
if (!dir.exists("logs"))   dir.create("logs")

# ---- Retrieve simulation index from command-line argument ----
args <- commandArgs(trailingOnly = TRUE)
idx  <- as.integer(args[1])

#------------------------------------------------------------------------------------------------------#
# ---- SHARED OUTPUT FILES (used across parallel processes) ----
#------------------------------------------------------------------------------------------------------#

# File paths for results and locks
lock_csv <- "./output/all_results.lock"
csv_file <- "./output/all_results.csv"

# Initialize the results CSV file if it does not already exist
lock <- lock(lock_csv)
if (!file.exists(csv_file)) {
  empty_dt <- data.table::data.table(
    sim                = integer(),
    year               = integer(),
    initial_population = numeric(),
    end_population     = numeric(),
    crossings          = numeric(),
    crossing_rate      = numeric(),
    collision_rate     = numeric(),
    death_collision    = numeric(),
    death_random       = numeric(),
    birth_total        = numeric(),
    net_growth         = numeric(),
    relative_growth    = numeric()
  )
  data.table::fwrite(empty_dt, csv_file)
}
unlock(lock)

# File paths for execution log
log_file <- "./output/simulation.txt"
lock_log <- "./output/simulation.lock"

# Log simulation start time
lock_l <- lock(lock_log)
write(paste0("Simulation ", idx, " started at ", Sys.time()),
      file = log_file, append = TRUE)
unlock(lock_l)

#------------------------------------------------------------------------------------------------------#
# ---- PARAMETER TABLE ----
# Defines all trait-by-traffic combinations passed to the IBM.
# 120 combinations: 4 categories x 3 crossing probability ranges x 4 traffic levels
#------------------------------------------------------------------------------------------------------#

# Road-crossing probability ranges (low / medium / high)
prob_cross <- dplyr::tibble(
  min_prob_cross = c(0,    0.33, 0.66),  # Lower bound of crossing probability
  max_prob_cross = c(0.33, 0.66, 1.00)   # Upper bound of crossing probability
)

# Build full parameter table using all combinations
pars <- tidyr::crossing(
  category     = 1:4,                    # Population categories (trait combinations)
  prob_cross,                            # Crossing probability ranges
  traffic_flow = c(0, 1, 5, 10)         # Traffic intensity levels (vehicles/min)
) |>
  dplyr::mutate(
    
    # ---- Movement Parameters ----
    steps     = 144,                                          # Number of movement steps per day
    mu        = 1,                                            # Shape parameter for truncated Pareto distribution
    daily_min = c(1000, 5000, 5000, 1000)[category],         # Minimum daily movement distance (m)
    daily_max = c(5000, 10000, 10000, 5000)[category],       # Maximum daily movement distance (m)
    lmin      = daily_min / steps,                           # Minimum step length (m)
    lmax      = daily_max / steps,                           # Maximum step length (m)
    angle_mean = 0,                                          # Mean turning angle (radians)
    angle_sd   = pi / 2,                                     # SD of turning angle (90 degrees of variability)
    
    # ---- Life History Parameters ----
    mature        = c(1, 1, 2, 3)[category],                 # Age at first reproduction (years)
    num_offspring = c(8, 3, 3, 2)[category],                 # Number of offspring per reproductive event
    reprod_period = 1,                                       # Reproductive period
    reprod_freq   = 1,                                       # Number of reproductive events per year
    K             = c(58, 55, 60, 95)[category],             # Carrying capacity
    
    # ---- Road Interaction Parameters ----
    speed_mean = 84,                                         # Mean individual speed when crossing the road (m/min)
    speed_sd   = 100,                                        # SD of individual speed when crossing the road (m/min)
    
    # Crossing probability bounds (inherited from prob_cross)
    min_prob_cross = min_prob_cross,                         # Minimum probability of attempting a road crossing
    max_prob_cross = max_prob_cross,                         # Maximum probability of attempting a road crossing
    
    # Traffic intensity (vehicles/min; used in collision_probability())
    traffic_flow = traffic_flow
    
  ) |>
  # Reorder columns for readability
  dplyr::relocate(min_prob_cross, .after = speed_sd)   |>
  dplyr::relocate(max_prob_cross, .after = min_prob_cross) |>
  dplyr::relocate(traffic_flow,   .after = max_prob_cross) |>
  dplyr::select(-daily_min, -daily_max)

# ---- Notes on parameter table ----
# 1. Categories represent theoretical populations with distinct trait combinations:
#      - Category-specific traits: step lengths (lmin/lmax), age at maturity,
#        and number of offspring per reproductive event.
# 2. traffic_flow modulates collision risk inside collision_probability().
# 3. Crossing probability bounds define individual behavioral tendency to
#    attempt road crossings (sampled uniformly between min and max per event).

#------------------------------------------------------------------------------------------------------#
# ---- PARALLEL SIMULATION EXECUTION ----
#------------------------------------------------------------------------------------------------------#

# Set up parallel processing using all available cores except one
future::plan(future::multisession, workers = parallel::detectCores() - 1)

# Start execution timer
tictoc::tic("Simulation time")

# Run IBM across all parameter combinations in parallel
go <- future.apply::future_lapply(
  1:nrow(pars),   # One iteration per parameter combination
  function(i) {
    
    # Log progress for the current iteration
    write(paste("Simulation", i, "/", nrow(pars), "- t:", Sys.time()),
          file = "./output/simulation.log", append = TRUE)
    
    # Extract parameter set for the current iteration
    row <- pars[i, ]
    
    # Run IBM with the selected parameter combination
    IBM(
      b          = 0.5,       # Base birth rate
      d          = 0.2,       # Base stochastic mortality rate
      K          = row$K,     # Carrying capacity
      N0         = 50,        # Initial population size
      days       = 365,       # Simulated days per year
      years      = 20,        # Total simulation duration (years)
      nrow       = 2000,      # Number of landscape rows (cells)
      ncol       = 2000,      # Number of landscape columns (cells)
      res        = 50,        # Cell resolution (meters)
      road_width = 5,         # Road width (meters)
      pars       = row,       # Full parameter set for this run
      sim        = idx,       # Simulation index (used as run identifier in output)
      lock_csv, csv_file, log_file, lock_log
    )
  },
  future.seed = 999           # Fixed seed for reproducibility across parallel runs
)

# Stop timer
tempo <- tictoc::toc()

#------------------------------------------------------------------------------------------------------#
# ---- FINALIZE AND SAVE RESULTS ----
#------------------------------------------------------------------------------------------------------#

# Log simulation completion time
lock_l <- lock(lock_log)
write(paste0("Simulation ", idx, " completed at ", Sys.time()),
      file = log_file, append = TRUE)
unlock(lock_l)

# Close parallel processing cluster
future::plan(future::sequential)

# Log save step
write(paste("Saving..."), file = "./output/simulation.log", append = TRUE)

# Consolidate results from all parallel runs into a single data.table
results_dt <- data.table::rbindlist(go, fill = TRUE)

# Save consolidated results as RDS (preserves R data types)
saveRDS(results_dt, file = "output/all_simulations_results.rds")

# Save consolidated results as CSV (for external use)
data.table::fwrite(results_dt, file = "output/all_simulations_results.csv")
