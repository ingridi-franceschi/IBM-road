#------------------------------------------------------------------------------------------------------#
#------------------------- INDIVIDUAL-BASED MODEL (IBM) -----------------------------------------------#
#------------------------------------------------------------------------------------------------------#
# Ecological traits and road traffic drive population persistence: a mechanistic modelling framework
#
# Date: 2026
# Version: 1.0
#------------------------------------------------------------------------------------------------------#
# DESCRIPTION:
#   This script implements an Individual-Based Model (IBM) that develops a 
#   mechanistic theoretical framework to explore how traffic flows influence 
#   population persistence across different ecological trait combinations.
#
#   The model pursues two main goals:
#     (1) To identify how populations differing in daily movement distance,
#         reproductive strategy, and road-crossing probability respond to
#         different traffic flows;
#     (2) To determine which trait combinations make a population more
#         vulnerable under different traffic scenarios.
#------------------------------------------------------------------------------------------------------#
# MODEL STRUCTURE:
#   Main model:
#     - IBM()                   : Main simulation loop (yearly > daily > steps)
#
#   Submodels:
#     - create_landscape()      : Builds the raster landscape
#     - create_individuals()    : Initializes individuals in habitat cells
#     - move_individuals()      : Lévy Walk movement in habitat
#     - collision_probability() : Estimates wildlife-vehicle collision risk
#     - random_death()          : Applies stochastic mortality
#     - birth()                 : Applies reproduction and offspring creation
#
#   Sub-submodels:
#     - move_habitat()          : Moves individual to adjacent habitat (avoidance)
#     - move_in_road()          : Moves individual onto the road (crossing attempt)
#     - move_out_road()         : Moves individual to the other side of the road (success crossing)
#------------------------------------------------------------------------------------------------------#
# INPUTS:
#   Population parameters : b (birth rate), d (death probability), K (carrying
#                           capacity), N0 (initial population size)
#   Time parameters       : days (per year), years (simulation duration)
#   Landscape parameters  : nrow, ncol (grid dimensions), res (cell resolution),
#                           road_width (meters)
#   Trait table           : pars (data frame with species traits and interaction
#                           parameters, including movement, speed, crossing
#                           probability, and reproduction)
#------------------------------------------------------------------------------------------------------#
# OUTPUTS:
# - Annual population size
# - Number of wildlife-vehicle collision and stochastic deaths
# - Birth number
# - Net growth and relative growth
#------------------------------------------------------------------------------------------------------#
# DEPENDENCIES:
#   terra        : Spatial raster operations and landscape management
#   data.table   : Fast data writing (fwrite) to CSV
#   truncnorm    : Truncated normal distribution for speed sampling
#   filelock     : File locking for parallel CSV writing
#
#   Install with: install.packages(c("terra", "data.table", "truncnorm", "filelock"))
#------------------------------------------------------------------------------------------------------#
# REFERENCE:
#   Litvaitis, J.A. & Tash, J.P. (2008). An approach toward understanding
#   wildlife-vehicle collisions. Environmental Management, 42, 688–697.
#   https://doi.org/10.1007/s00267-008-9108-4
#------------------------------------------------------------------------------------------------------#

# Clean environment
rm(list=ls())

########### MODEL IBM ########### 
IBM <- function(# ---- Population Parameters ----
                b,          # Base birth rate
                d,          # Base death probability
                K,          # Carrying capacity
                N0,         # Initial population size  
                # ---- Time Parameters ----
                days,       # Number of days per year
                years,      # Number of simulation years
                # ---- Landscape Parameters ----
                nrow, ncol, # Landscape dimensions (in cells)
                res,        # Cell resolution
                road_width,   # Road width (meters)
                # ---- Parameters ----
                pars,        # Table containing population traits and environmental interactions
                # ---- Simulation Information ----
                sim
) {
  
  # ---- Initialize Landscape ----
  # Creates a raster with habitat (value = 1), risk areas (value = 2), and roads (value = 3)
  env <- create_landscape(nrow = nrow, ncol = ncol, res = res, pars = pars, plot = FALSE)
  
  # ---- Initialize Individuals ----
  # Creates N0 individuals with random positions in habitat
  indiv <- create_individuals(N0 = N0, env = env, pars = pars, plot_point = FALSE)
  
  # ---- Data Structures for Results ----
  # Track dead individuals (separated by mortality cause)
  indiv_dead <- list(collision = list(), random = list())
  
  # Vectors to store yearly data
  initial_population <- numeric(years)  # Vector to store the initial population size for each year
  end_population     <- numeric(years)  # Vector to store the final population size for each year
  death_collision    <- numeric(years)  # Vector to store the number of wildlife-vehicle collision per year
  death_random       <- numeric(years)  # Vector to store the number of stochastic deaths per year
  birth_total        <- numeric(years)  # Vector to store the total births per year
  crossings          <- numeric(years)  # Vector to store the number of crossings per year
  crossing_rate      <- numeric(years)  # Vector to store the crossing rate per year
  collision_rate     <- numeric(years)  # Vector to store the wildlife-vehicle collision rate per year
  net_growth         <- numeric(years)  # Vector to store net population growth (births - deaths) per year
  relative_growth    <- numeric(years)  # Vector to store relative growth rate (net_growth / initial_population) per year
  collision_start_year <- 6 # Year that starts vehicle collision
  
  # ---- Yearly Loop ----
  for (year in 1:years) { 
    if(length(indiv) > 0) {
      # Store initial population size at the start of the year (counts alive individuals)
      initial_population[year] <- sum(vapply(indiv, function(x) x$alive == 1, logical(1)))
      
      indiv_dead_ids <- c() # Vector to track IDs of individuals dead by collision daily
      cross_count <- numeric(length(indiv)) # Vector to count crossings per individual
      
      # ---- Daily Loop ----
      for (day in 1:days) { 
        # --- write log global ---
        lock_l <- lock(lock_log)
        write(paste0("Simulation",sim,"- Year: ",year,"- Day: ",day,"- t:",Sys.time()), 
              file = log_file, append = TRUE)
        unlock(lock_l)
        # ---- Individual Movement ----
        for (i in seq_along(indiv)) {
          ind <- indiv[[i]]
          
          # ---- Movement Steps Loop ----
          for (step in seq_len(pars$steps)) { 
            # Get current cell type (1 = habitat, 2 = risk, 3 = road)
            i_cell <- terra::cellFromXY(env, cbind(ind$x, ind$y))
            cell_type <- env[i_cell]
            
            # Case 1: Individual is in a habitat cell
            if (cell_type == 1) { 
              ind <- move_individuals(env, ind, pars, plot_point = FALSE)
            } 
            # Case 2: Individual is in a risk area - decide to cross the road or avoid
            else if (cell_type == 2) {
              # Calculate the probability of crossing
              cross_prob <- runif(1, pars$min_prob_cross, pars$max_prob_cross)
              # Decision: Cross the road (move to road cell)
              if (cross_prob > runif(1)) {
                ind <- move_in_road(env, ind, plot_point = FALSE) # Cross into road
                cross_count[[i]] <- cross_count[[i]] + 1 # Increment crossing count
              } else  {
                # Decision: Avoid crossing (move to safe habitat cell)
                ind <- move_habitat(env, ind, pars, plot_point = FALSE)
              }
            }
            
            # Case 3: Individual in a road cell - risk of wildlife-vehicle collision
            else if(cell_type == 3){
              if (year >= collision_start_year) {
                # Calculate the probability of wildlife-vehicle collision
                collision_result <- collision_probability(env, ind, pars, road_width)
                # Check if individual died from vehicle collision
                if (collision_result$alive == 0 && collision_result$death_type == "collision") {
                  # Ensure the individual isn't already recorded as dead (avoid duplicates)
                  if (!(ind$ID %in% indiv_dead_ids)) {
                    # Add dead individual to the roadkill death list
                    indiv_dead$collision <- append(indiv_dead$collision, list(collision_result))
                    # Track the ID to prevent duplicate records
                    indiv_dead_ids <- c(indiv_dead_ids, ind$ID)
                  }
                  break  # Stop movement if dead
                } else {
                  # If survives, move out of the road
                  ind <- move_out_road(env, res, ind, pars, plot_point = FALSE)
                }
              } else
                ind <- move_out_road(env, res, ind, pars, plot_point = FALSE)
            }
          } # End steps loop
          indiv[[i]] <- ind  # Update individual after all movements
        } # End individuals loop
        
        rm(ind) # Free memory
        
        # Remove of roadkill individuals from the population
        dead_collision_ids <- unique(unlist(lapply(indiv_dead$collision, function(y) y$ID)))
        indiv <- Filter(function(x) !(x$ID %in% dead_collision_ids), indiv)
        
        # Reset daily movement tracks for all survivors
        indiv <- lapply(indiv, function(a) { a$track <- a$track[nrow(a$track) -3:0, 1:2]; return(a)})
        
      } # End daily loop
      
      if(length(indiv) > 0) {
        # ---- Stochastic Mortality ----
        # Apply random death probability to all individuals alive
        result_death_random <- random_death(indiv, d)
        # Add stochastically dead individuals to the random mortality list
        indiv_dead$random <- c(indiv_dead$random, Filter(function(x) x$alive == 0, result_death_random$indiv))
        # Record stochastic death count for the year
        death_random[year] <- result_death_random$d_random_count 
        
        # Remove dead individuals from the population
        dead_random_ids <- unique(unlist(lapply(indiv_dead$random, function(y) y$ID)))
        indiv <- Filter(function(x) !(x$ID %in% dead_random_ids), indiv)
        
        # ---- Reproduction ----
        # Calculate new births
        result_birth <- birth(env, indiv, indiv_dead, pars, b, K)
        # Store newborn individuals
        new_indiv <- result_birth$indiv 
        
        # ---- Population Updates ----
        # Age individuals
        indiv <- lapply(indiv, function(x) { x$age <- x$age + 1; return(x) })
        # Add newborns to the population
        indiv <- c(indiv, new_indiv)
        
        # ---- Annual Metrics ----

        # Count wildlife-vehicle collision for the year
        death_collision[year] <- length(indiv_dead_ids)
        # Wildlife-vehicle collision rate
        collision_rate[year] <- death_collision[year] / days
        # Count total road crossings for the year
        crossings[year] <- sum(cross_count)
        # Crossing rate
        crossing_rate[year] <- crossings[year] / days
        # Record total births for the year
        birth_total[year] <- result_birth$birth_count
        # Final population size
        end_population[year] <- initial_population[year] - death_random[year] - death_collision[year] + birth_total[year]
        # Net population change
        net_growth[year] <- end_population[year] - initial_population[year]
        # Relative growth rate
        relative_growth[year] <- net_growth[year] / initial_population[year]
      }
      
      # ---- Save annual result immediately ----
      result_year <- data.table::data.table(
        sim = sim,
        year = year,
        initial_population = initial_population[year],
        end_population = end_population[year],
        crossings = crossings[year],
        crossing_rate = crossing_rate[year],
        collision_rate = collision_rate[year],
        death_collision = death_collision[year],
        death_random = death_random[year],
        birth_total = birth_total[year],
        net_growth = net_growth[year],
        relative_growth = relative_growth[year]
      )
      
      lock <- lock(lock_csv)
      data.table::fwrite(result_year, csv_file, append = TRUE)
      unlock(lock)
    }  
    else{
      break # Terminate simulation if population goes extinct
    }  
  } # End yearly loop
  gc() # Clean memory
} # End IBM model

################################################################################
# ********************************* SUBMODELS *********************************#
################################################################################
# Function to create the landscape (grid) with habitat, road, and risk area
create_landscape <- function(nrow, ncol, res, pars, plot = TRUE) {
  
  # Create a raster object with specified number of rows (nrow), columns (ncol), and resolution (res)
  env <- terra::rast(nrows = nrow, ncols = ncol, 
                     xmin = 0, xmax = ncol * res, 
                     ymin = 0, ymax = nrow * res)
  
  # Initialize all raster cells with value 1, representing habitat
  env[] <- 1
  
  # Define the road in the middle of the landscape (nrow/2) and in all columns (1:ncol), and define value 3
  road_area <- nrow %/% 2
  road_cells <- ((road_area - 1) * ncol + 1):(road_area * ncol)
  env[road_cells] <- 3
  
  # Select the maximium step length (lmax) of the population to define the risk area
  margin_width <- ceiling(pars$lmax/res)
  
  # Define the risk area cells
  risk_area <- road_cells
  for (i in 1:margin_width) {
    risk_area <- c(risk_area, road_cells - i * ncol(env), road_cells + i * ncol(env))
  }
  
  # Define the road margins as risk areas with value 2
  env[risk_area[which(env[risk_area] == 1)]] <- 2
  
  # Set colors for habitat (dark green), margins (olive green), and road (gray)
  cores <- c("darkgreen", "darkolivegreen", "darkgray")
  
  # Plot the landscape with defined colors
  if(plot == TRUE){
    terra::plot(env, col = cores)
  }
  # Return the raster representing the landscape (with habitat, risk area, and road)
  return(env)
}

# Function to create initial individuals in the landscape
create_individuals <- function(N0, env, pars, plot_point = TRUE) {
  # Initialize empty list to store individuals (N0)
  indiv <- vector(mode = "list", N0)
  
  # Identify cells classified as habitat (value == 1)
  habitat_cells <- which(env[] == 1)
  
  # Loop to create each individual
  for (i in seq(indiv)) {
    cell <- sample(habitat_cells, 1) # Randomly select a habitat cell
    coords <- terra::xyFromCell(env, cell) # Get coordinates of the selected cell
    
    # Individual attributes
    indiv[[i]]$ID <- i # Unique ID
    indiv[[i]]$alive <- 1  # Individual starts alive
    indiv[[i]]$age <- sample(0:3, 1, replace = TRUE)  # Random initial age (0 = juvenile, 1–2 = adult)
    indiv[[i]]$angle_abs <- rnorm(1, mean = pars$angle_mean, sd = pars$angle_sd) # Initial angle
    indiv[[i]]$color <- "blue"  # Initial color
    
    # Define x and y coordinates to the individual based on the selected cell
    indiv[[i]]$x <- coords[1]
    indiv[[i]]$y <- coords[2]
    
    # Create a data frame to track individual movement during daily loops
    indiv[[i]]$track <- data.frame(x = indiv[[i]]$x, y = indiv[[i]]$y)
    
    # Optionally plot the individual
    if(plot_point == TRUE){
      points(indiv[[i]]$x, indiv[[i]]$y, col = indiv[[i]]$color, pch = 19, cex = 0.5)
    }
  }
  # Return the list of individuals
  return(indiv)
}

# Function to move individuals using the Lévy Walk based on truncaded Pareto distribution with a toroidal landscape
move_individuals <- function(env, indiv, pars, plot_point = TRUE) {
  # Determine step length using truncated Pareto distribution
  step_length <- rpareto_trunc(1, pars$lmin, pars$lmax, pars$mu)
  
  # Get individual's current angle
  current_angle <- indiv$angle_abs
  
  # Calculate a new angle based on normal distribution, it defines the movement direction
  new_angle <- rnorm(1, mean = pars$angle_mean, sd = pars$angle_sd) %% (2 * pi)
  
  # Calculate a turning angle considering the current and the new angle
  turning_angle <- (current_angle + new_angle) %% (2 * pi)
  
  # Update individual's current angle
  indiv$angle_abs <- turning_angle
  
  # Calculate new X and Y positions based on angle and step length
  new_x <- indiv$x + cos(turning_angle) * step_length
  new_y <- indiv$y + sin(turning_angle) * step_length
  
  # Apply toroidal boundaries for x and y positions
  new_x <- ((new_x - terra::xmin(env)) %% terra::xmax(env)) + terra::xmin(env)
  new_y <- ((new_y - terra::ymin(env)) %% terra::ymax(env)) + terra::ymin(env)
  
  # Update individual's position
  indiv$x <- new_x
  indiv$y <- new_y
  
  # Add new position to individual's track
  indiv$track <- rbind(indiv$track, 
                       data.frame(x = new_x, y = new_y))
  
  # Optionally plot the individual's position after moving
  if(plot_point == TRUE) {
    points(indiv$x, indiv$y, col = "orange", pch = 19, cex = 0.5)
  }
  # Return the updated individual
  return(indiv)
}

# Function to evaluate the probability of wildlife-vehicle collision
collision_probability <- function(env, indiv, pars, road_width) {
  
  # Calculate the individual's speed based on a truncated normal distribution with defined minimum and maximum speed
  speed <- truncnorm::rtruncnorm(n = 1, a = 50, b = 700, mean = pars$speed_mean, sd = pars$speed_sd)
  
  # Calculate the time the individual takes to cross the road
  time_cross <- road_width / speed
  
  # Estimate collision probability based on Litvaitis and Tash (2008)
  collision_prob <- 1 - exp(-pars$traffic_flow * time_cross)
  
  # Determine if a collision occurs (individual dies if collision probability is higher than a random draw)
  if (collision_prob > runif(1)) {
    indiv$alive <- 0  # Set alive status to 0
    indiv$collision <- collision_prob # Store the collision probability value
    indiv$death_type <- "collision"  # Set the death type as collision
  } 
  # Return the updated individual  
  return(indiv)
}

# Function to apply stochastic mortality
random_death <- function(indiv, d) {
  d_random_count <- 0 # Counter for stochastic death
  
  # For an alive individual
  for (i in seq_along(indiv)) {
    # Apply stochastic death based on probability d
    if (indiv[[i]]$alive == 1 && rbinom(1, 1, d) == 1) { 
      # If the result is 1, the individual dies
      indiv[[i]]$alive <- 0  # Set alive status to 0
      indiv[[i]]$death_type <- "random"  # Set the death type as random
      d_random_count <- d_random_count + 1 # Increment the stochastic death counter
    }
  } 
  # Return a list containing updated individuals and stochastic death count
  return(list(indiv = indiv, d_random_count = d_random_count))
}

# Function for births
birth <- function(env, indiv, indiv_dead, pars, b, K) {
  birth_count <- 0 # Counter for number of births
  new_indiv_list <- list() # List to store newborn individuals 
  
  # Combine living and dead individuals to maintain global ID tracking
  indiv_total = c(indiv, indiv_dead$collision, indiv_dead$random)
  # Get the highest existing ID among all individuals (living and dead) to ensure unique IDs
  ID_number <- max(sapply(indiv_total, function(x) x$ID))
  
  # Iterate over all living individuals to calculate birth probability
  for (j in seq_along(indiv)) { 
    # Check if the individual has reached reproductive age and if it's a reproductive period
    if(indiv[[j]]$age >= pars$mature && (indiv[[j]]$age %% pars$reprod_period == 0)) {
      # Calculate adjusted birth probability based on reproductive frequency
      new_birth_prob <- 1 - (1 - b)^pars$reprod_freq
      # Adjust the probability based on population density and carrying capacity
      new_b <- max(0, min(new_birth_prob * (1 - length(indiv)/K), 1)) # Limit probability between 0 and 1
      
      # Create the total number of offspring based on Poisson distribution
      offspring <- rpois(1, lambda = new_b * pars$num_offspring)
      
      # If one or more new individuals are created
      if (offspring > 0) {
        # Create new individuals based on the number of offspring
        for (f in 1:offspring) {  
          habitat_cells <- which(env[] == 1) # Select random habitat cells in the landscape (value = 1)
          cell <- sample(habitat_cells, 1) # Randomly choose a habitat cell
          coords <- terra::xyFromCell(env, cell) # Get X and Y coordinates of the selected cell
          ID_number <- ID_number + 1 # Update the new individual's ID
          
          # Create the new individual with attributes
          new_indiv <- list(
            ID = ID_number, # Unique ID
            alive = 1, # Individual starts alive
            age = 0, # Initial age of the newborn is 0
            angle_abs = rnorm(1, mean = pars$angle_mean, sd = pars$angle_sd), # Initial movement angle
            x = coords[1], # X coordinate
            y = coords[2], # Y coordinate
            track = data.frame(x = coords[1], y = coords[2]) # Movement track
          ) 
          # Add the new individual to the list
          new_indiv_list <- c(new_indiv_list, list(new_indiv))
          # Increment birth counter
          birth_count <- birth_count + 1  
        }
      }
    }
  }  
  # Return the list of new individuals and total births
  return(list(indiv = new_indiv_list, birth_count = birth_count))  
}

################################################################################
# ******************************* SUB-SUBMODELS *******************************#
################################################################################
# Function to move an individual to an adjacent habitat cell. This is used when the individual decides not to cross the road
move_habitat <- function(env, indiv, pars, plot_point = TRUE) {
  
  # Get the coordinates of the individual's previous cell (last point in its track)
  previous_cel <- indiv$track[nrow(indiv$track) -1, 1:2]
  
  # Get the cell index for the last known coordinates
  last_cell <- terra::cellFromXY(env, cbind(previous_cel$x, previous_cel$y))
  
  # Identify all 8 adjacent neighbor cells
  neig <- as.vector(terra::adjacent(env, last_cell, directions = 8))
  
  # Filter for valid habitat cells (value == 1) and randomly select one
  new_cell <- sample(neig[!is.na(terra::values(env)[neig]) & terra::values(env)[neig] == 1], 1)
  
  # Get x and y coordinates of the selected cell
  new_coords <- terra::xyFromCell(env, new_cell)
  
  # Define x and y coordinates to the individual based on the selected cell
  indiv$x <- new_coords[1]
  indiv$y <- new_coords[2]
  
  # Add new position to individual's track
  indiv$track <- rbind(indiv$track, data.frame(x = indiv$x, y = indiv$y))
  
  # Optionally plot the individual's position after moving
  if(plot_point == TRUE){
    points(indiv$x, indiv$y, col = "red", pch = 19, cex = 0.5)
  }
  # Return the updated individual
  return(indiv)
}

# Function to move individuals into road after deciding to cross the road
move_in_road <- function(env, indiv, plot_point = TRUE){
  
  # Get the index of the current cell from x and y coordinates
  current_xy <- cbind(indiv$x, indiv$y)
  
  # Define a new position on the road
  road_y <- terra::yFromRow(env, (nrow(env) / 2))
  
  # Update individual's position
  indiv$x <- current_xy[1] # Keep the same column
  indiv$y <- road_y # Position in the road row
  
  # Add new position to individual's track
  indiv$track <- rbind(indiv$track, data.frame(x = indiv$x, y = indiv$y))
  
  # Optionally plot the individual's position after moving
  if(plot_point == TRUE){
    points(indiv$x, indiv$y, col = "black", pch = 19, cex = 0.5)
  }
  # Return the updated individual
  return(indiv)
}

# Function for individuals to cross the road
move_out_road <- function(env, res, indiv, pars, plot_point = TRUE) {
  # Get the index of the current cell from x and y coordinates
  current_cell <- terra::cellFromXY(env, cbind(indiv$x, indiv$y))
  
  # Get the coordinates of the individual's previous cell (last point in its track)
  previous_xy <- indiv$track[nrow(indiv$track) -1, 1:2]
  
  # Get the index of the previous cell from track coordinates
  previous_cel <- terra::cellFromXY(env, cbind(previous_xy$x, previous_xy$y))
  
  # Select the maximium step length (lmax) of the population to define the risk area
  margin_width <- ceiling(pars$lmax/res)
  
  # Define the road position
  mid_row <- nrow(env) / 2
  
  # Determine if the last individual position was located above or below the road
  current_row <- terra::rowColFromCell(env, previous_cel)[1]
  
  if(current_row < mid_row){
    # If the individual was above the road, move it to cross above the road in a habitat cell
    new_cell <- current_cell + ncol(env) * (margin_width + 1)
  } else { 
    # Otherwise, move the individual to cross below the road in a habitat cell
    new_cell <- current_cell - ncol(env) * (margin_width + 1)
  }
  
  # Get x and y coordinates of the selected cell
  new_coords <- terra::xyFromCell(env, new_cell)
  
  # Get individual's current angle
  current_angle <- indiv$angle_abs
  
  # Calculate a new angle based on normal distribution, it defines the movement direction
  new_angle <- (rnorm(1, mean = pars$angle_mean, sd = pars$angle_sd)) %% (2 * pi) 
  
  # Calculate a turning angle considering the current and the new angle
  turning_angle <- (current_angle + new_angle + pi) %% (2 * pi)
  
  # Update individual's current angle
  indiv$angle_abs <- turning_angle
  
  # Calculate new X and Y positions based on turning angle and step length
  new_x <- new_coords[1] + cos(turning_angle)
  new_y <- new_coords[2] + sin(turning_angle)
  
  # Apply toroidal boundaries for x and y positions
  new_x <- ((new_x - terra::xmin(env)) %% terra::xmax(env)) + terra::xmin(env)
  new_y <- ((new_y - terra::ymin(env)) %% terra::ymax(env)) + terra::ymin(env)
  
  # Update individual's position
  indiv$x <- new_x
  indiv$y <- new_y
  
  # Add new position to individual's track
  indiv$track <- rbind(indiv$track, data.frame(x = new_x, y = new_y))
  
  # Optionally plot the individual's position after moving
  if(plot_point == TRUE){
    points(indiv$x, indiv$y, col = "purple", pch = 19, cex = 0.5)
  }
  # Return the updated individual
  return(indiv)
}

######## TRUNCATED PARETO DISTRIBUTION FUNCTION ########
rpareto_trunc <- function(num_step, lmin, lmax, mu) {
  u <- runif(num_step)
  ((lmax^(-mu) - lmin^(-mu)) * u + lmin^(-mu))^(-1 / mu)
}