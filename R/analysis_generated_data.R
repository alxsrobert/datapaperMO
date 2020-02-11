## Generate all o2geosocial runs on simulated dataset (toy_outbreak_long)

source("R/library_importation.R")
source("R/function_generate_dataset.R")
source("R/load_data_distributions.R")


#### Define parameters and generate dataset, or use default toy_outbreak_long dataset ####

# If FALSE: Use the data attached to o2geosocial, if TRUE: generate the data
generate <- TRUE

if(generate){
  # Set seed
  set.seed(1)
  
  # Draw the r0 in each state, following an exponential distribution
  r0_state <- rexp(n = length(unique(dt_state_county$STNAME)), 
                   rate = 1.5)
  names(r0_state) <- unique(dt_state_county$STNAME)
  r0_state[r0_state>1] <- 0.9
  
  # Spatial parameters
  c <- 0.84
  b <- 0.14
  
  # Spatial threshold
  gamma <- 100
  
  # Maximum number of cases in the generated dataset
  nb_cases <- 2200
  
  toy_outbreak_long <- generate_dataset(b = b, c = c, gamma = gamma, 
                                        dt_distance = dt_distance,
                                        polymod_prop = polymod_prop, w = w, 
                                        nb_cases = nb_cases, r0_state = r0_state, 
                                        pop_county = pop_county, 
                                        dt_state_county = dt_state_county,
                                        t_min = as.Date("2003-01-01"), 
                                        t_max = as.Date("2016-01-01"))
  states_keep <- names(table(toy_outbreak_long[["cases"]]$State)[table(toy_outbreak_long[["cases"]]$State) > 20])
  toy_outbreak_long[["cases"]] <- toy_outbreak_long[["cases"]][is.element(State, states_keep),]
} else
  data("toy_outbreak_long")

# Save elements of toy_outbreak_long
dt_cases <- toy_outbreak_long[["cases"]]
dt_cases <- dt_cases[order(Date), ]
dt_regions <- toy_outbreak_long[["dt_regions"]]
all_dist <- distGeo(matrix(c(rep(dt_regions$long, nrow(dt_regions)), 
                             rep(dt_regions$lat, nrow(dt_regions))), ncol = 2), 
                    matrix(c(rep(dt_regions$long, each = nrow(dt_regions)), 
                             rep(dt_regions$lat, each = nrow(dt_regions))), ncol = 2))

dist_mat <- matrix(all_dist/1000, nrow = nrow(dt_regions))
pop_vect <- dt_regions$population
age_contact <- toy_outbreak_long[["age_contact"]]

names(pop_vect) <- rownames(dist_mat) <- colnames(dist_mat) <- dt_regions$region

f_null <- function(data, param) {
  return(0.0)
}
#### o2geosocial analysis toy data: threshold and import ####

## No import, threshold = 0.01
data <- outbreaker_data(dates = dt_cases$Date,
                        age_group = dt_cases$age_group,
                        region = dt_cases$county,
                        w_dens = dnorm(x = 1:300, mean = 11.7, sd = 2.0),
                        f_dens = dgamma(x = 1:300, scale = 0.4286, shape = 26.834),
                        a_dens = as.matrix(age_contact),
                        genotype = dt_cases$Genotype,
                        is_cluster = NULL,
                        population = pop_vect,
                        distance = dist_mat,
                        import = FALSE)

config <- create_config(data = data, 
                        init_tree = "star", 
                        spatial_method = "exponential",
                        gamma = 100,
                        delta = 30,
                        n_iter = 50000, 
                        sample_every = 50, 
                        max_kappa = 5,
                        find_import = TRUE, 
                        outlier_threshold = 0.01, 
                        outlier_relative = FALSE,
                        n_iter_import = 15000, 
                        sample_every_import = 50, verbatim = TRUE,
                        burnin = 5000)
priors <- custom_priors()
likelihoods <- custom_likelihoods()
moves <- custom_moves()
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)

saveRDS(out, file = "toy_outbreak_runs/no_import_thresh001.rds")

## No import, threshold = 0.05
config$outlier_threshold <- 0.05
config <- create_config(config, data = data)
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
saveRDS(out, file = "toy_outbreak_runs/no_import_thresh005.rds")

## No import, threshold = 0.9 (Relative)
config$outlier_threshold <- 0.9
config$outlier_relative <- TRUE
config <- create_config(data = data, 
                        config = config)


out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)

saveRDS(out, file = "toy_outbreak_runs/no_import_thresh09.rds")

## No import, threshold = 0.95 (Relative)
config$outlier_threshold <- 0.95
config <- create_config(data = data, 
                        config = config)


out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)

saveRDS(out, file = "toy_outbreak_runs/no_import_thresh095.rds")

## Import, no further estimation of import status
data$import <- dt_cases$import
config$find_import <- FALSE
data <- outbreaker_data(data = data)
config <- create_config(data = data, 
                        config = config)
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
saveRDS(out, file = "toy_outbreak_runs/with_import.rds")

## Import + threshold = 0.05
config$find_import <- TRUE
config$outlier_threshold <- 0.05
config$outlier_relative <- FALSE
config <- create_config(data = data, 
                        config = config)
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
saveRDS(out, file = "toy_outbreak_runs/with_import_005.rds")

## Import + threshold = 0.95 (Relative)
data <- outbreaker_data(data = data,
                        import = dt_cases$import)
config$outlier_threshold <- 0.95
config$outlier_relative <- TRUE
config <- create_config(data = data, 
                        config = config)


out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
saveRDS(out, file = "toy_outbreak_runs/with_import_095.rds")
stop()
#### o2geosocial analysis toy data: elements of likelihood ####

## Import, time only
data <- outbreaker_data(dates = dt_cases$Date,
                        age_group = dt_cases$age_group,
                        region = dt_cases$county,
                        f_dens = matrix(exp(f), nrow = 1),
                        w_dens = matrix(exp(w), nrow = 1), 
                        a_dens = NULL, 
                        genotype = rep("Not attributed", dim(dt_cases)[1]),
                        is_cluster = NULL,
                        population = pop_vect,
                        distance = dist_mat,
                        import = dt_cases$import)

config <- create_config(data = data, max_kappa = 2,
                        init_tree = "star", 
                        spatial_method = "exponential",
                        init_kappa = 1,
                        gamma = NULL,
                        delta = 30,
                        n_iter = 50000, 
                        n_iter_import = 15000, 
                        sample_every = 50, 
                        sample_every_import = 50, 
                        burnin = 5000,
                        max_kappa = 5,
                        find_import = FALSE,
                        outlier_threshold = 0.05, verbatim = T,
                        outlier_relative = FALSE,
                        move_a = FALSE,
                        move_b = FALSE,
                        move_pi = FALSE,
                        move_t_inf = TRUE,
                        move_kappa = TRUE,
                        move_swap_cases = TRUE)
priors <- custom_priors()
likelihoods <- custom_likelihoods(space = f_null, 
                                  age = f_null)
moves <- custom_moves()
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)

saveRDS(out, file = paste0("toy_outbreak_runs/import_time_only.rds"))

## Import, time and genotype
data$genotype <- dt_cases$Genotype
data <- outbreaker_data(data = data)
config$init_alpha <- NULL
config <- create_config(data = data,
                        config)
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
saveRDS(out, file = paste0("toy_outbreak_runs/import_time_genotype.rds"))


## Import, time, genotype and space
likelihoods <- custom_likelihoods(age = f_null)
config$move_a <- TRUE
config$move_b <- TRUE
config$gamma <- 100
config$init_alpha <- NULL

config <- create_config(data = data,
                        config)
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
saveRDS(out, file = paste0("toy_outbreak_runs/import_time_genotype_space.rds"))


## Import, time, genotype and age
data$a_dens = matrix(polymod_prop, nrow = nrow(polymod_prop))
data <- outbreaker_data(data = data)
config$move_a <- FALSE
config$move_b <- FALSE
config$gamma <- NULL
config <- create_config(data = data,
                        config)
likelihoods <- custom_likelihoods(space = f_null)
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
saveRDS(out, file = paste0("toy_outbreak_runs/import_time_genotype_age.rds"))


#### o2geosocial analysis toy data: Different genotype report rate ####

# Set seed
set.seed(1)

# Draw the r0 in each state, following an exponential distribution
r0_state <- rexp(n = length(unique(dt_state_county$STNAME)), 
                 rate = 1.5)
names(r0_state) <- unique(dt_state_county$STNAME)
r0_state[r0_state>1] <- 0.9

# Spatial parameters
c <- 0.84
b <- 0.14

# Spatial threshold
gamma <- 100

# Maximum number of cases in the generated dataset
nb_cases <- 2200

toy_outbreak_long <- generate_dataset(b = b, c = c, gamma = gamma, 
                                      dt_distance = dt_distance,
                                      polymod_prop = polymod_prop, w = w, 
                                      nb_cases = nb_cases, r0_state = r0_state, 
                                      pop_county = pop_county, 
                                      dt_state_county = dt_state_county,
                                      t_min = as.Date("2003-01-01"), 
                                      t_max = as.Date("2016-01-01"))
states_keep <- names(table(toy_outbreak_long[["cases"]]$State)[table(toy_outbreak_long[["cases"]]$State) > 20])
toy_outbreak_long[["cases"]] <- toy_outbreak_long[["cases"]][is.element(State, states_keep),]

## No import, threshold = 0.05, all genotype reported
data <- outbreaker_data(dates = dt_cases$Date,
                        age_group = dt_cases$age_group,
                        region = dt_cases$county,
                        w_dens = dnorm(x = 1:300, mean = 11.7, sd = 2.0),
                        f_dens = dgamma(x = 1:300, scale = 0.4286, shape = 26.834),
                        a_dens = as.matrix(age_contact),
                        genotype = dt_cases$Genotype,
                        is_cluster = NULL,
                        population = pop_vect,
                        distance = dist_mat,
                        import = FALSE)

config <- create_config(data = data, 
                        init_tree = "star", 
                        spatial_method = "exponential",
                        gamma = 100,
                        delta = 30,
                        n_iter = 50000, 
                        sample_every = 50, 
                        max_kappa = 5,
                        find_import = TRUE, 
                        outlier_threshold = 0.05, 
                        outlier_relative = FALSE,
                        n_iter_import = 15000, 
                        sample_every_import = 50, 
                        burnin = 5000)
priors <- custom_priors()
likelihoods <- custom_likelihoods()
moves <- custom_moves()

out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)

saveRDS(out, file = "toy_outbreak_runs/no_import_thresh005_all_gen.rds")


# Set seed
set.seed(1)

# Draw the r0 in each state, following an exponential distribution
r0_state <- rexp(n = length(unique(dt_state_county$STNAME)), 
                 rate = 1.5)
names(r0_state) <- unique(dt_state_county$STNAME)
r0_state[r0_state>1] <- 0.9

# Spatial parameters
c <- 0.84
b <- 0.14

# Spatial threshold
gamma <- 100

# Maximum number of cases in the generated dataset
nb_cases <- 2200

toy_outbreak_long <- generate_dataset(b = b, c = c, gamma = gamma, 
                                      dt_distance = dt_distance,
                                      polymod_prop = polymod_prop, w = w, 
                                      nb_cases = nb_cases, r0_state = r0_state, 
                                      pop_county = pop_county, 
                                      dt_state_county = dt_state_county,
                                      t_min = as.Date("2003-01-01"), 
                                      t_max = as.Date("2016-01-01"))
states_keep <- names(table(toy_outbreak_long[["cases"]]$State)[table(toy_outbreak_long[["cases"]]$State) > 20])
toy_outbreak_long[["cases"]] <- toy_outbreak_long[["cases"]][is.element(State, states_keep),]

## No import, threshold = 0.05, no genotype reported
data <- outbreaker_data(dates = dt_cases$Date,
                        age_group = dt_cases$age_group,
                        region = dt_cases$county,
                        w_dens = dnorm(x = 1:300, mean = 11.7, sd = 2.0),
                        f_dens = dgamma(x = 1:300, scale = 0.4286, shape = 26.834),
                        a_dens = as.matrix(age_contact),
                        genotype = dt_cases$Genotype,
                        is_cluster = NULL,
                        population = pop_vect,
                        distance = dist_mat,
                        import = FALSE)

config <- create_config(data = data, 
                        init_tree = "star", 
                        spatial_method = "exponential",
                        gamma = 100,
                        delta = 30,
                        n_iter = 50000, 
                        sample_every = 50, 
                        max_kappa = 5,
                        find_import = TRUE, 
                        outlier_threshold = 0.05, 
                        outlier_relative = FALSE,
                        n_iter_import = 15000, 
                        sample_every_import = 50, 
                        burnin = 5000)
priors <- custom_priors()
likelihoods <- custom_likelihoods()
moves <- custom_moves()

out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)

saveRDS(out, file = "toy_outbreak_runs/no_import_thresh005_no_gen.rds")
