## Generate all measlesoutbreaker runs on simulated dataset (toy_outbreak)

source("R/library_importation.R")
source("R/function_generate_dataset.R")
source("R/load_data_distributions.R")


#### Define parameters and generate dataset, or use default toy_outbreak dataset ####

# If FALSE: Use the data attached to measlesoutbreaker, if TRUE: generate the data
generate <- TRUE

if(generate){
  # Set seed
  set.seed(1)
  
  # Draw the r0 in each state, following an exponential distribution
  r0_state <- rexp(n = length(unique(dt_state_county$STNAME)), 
                   rate = 1.5)
  names(r0_state) <- unique(dt_state_county$STNAME)
  r0_state[r0_state>1] <- 0.95
  
  # Spatial parameters
  a <- 0.84
  b <- 0.08
  
  # Spatial threshold
  gamma <- 100
  
  # Maximum number of cases in the generated dataset
  nb_cases <- 1000
  
  toy_outbreak <- generate_dataset(a = a, b = b, gamma = gamma, 
                                   dt_distance = dt_distance,
                                   polymod_prop = polymod_prop, w = w, 
                                   nb_cases = nb_cases, r0_state = r0_state, 
                                   pop_county = pop_county, 
                                   dt_state_county = dt_state_county,
                                   t_min = as.Date("2010-01-01"), 
                                   t_max = as.Date("2016-01-01")) 
} else
  data("toy_outbreak")
# Save elements of toy_outbreak
dt_cases <- toy_outbreak[["cases"]]
dt_cases <- dt_cases[order(Date), ]
dist_mat <- toy_outbreak[["distance"]]
pop_vect <- toy_outbreak[["population"]]
age_contact <- toy_outbreak[["age_contact"]]


f_null <- function(data, param) {
  return(0.0)
}


#### measlesoutbreaker analysis toy data: threshold and import ####

## No import, threshold = 0.01
data <- outbreaker_data(dates = dt_cases$Date,
                        age_group = dt_cases$age_group,
                        postcode = dt_cases$county,
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
                        n_iter = 70000, 
                        sample_every = 50, 
                        min_date = -10, 
                        max_kappa = 5,
                        find_import = TRUE, 
                        outlier_threshold = 0.01, 
                        outlier_relative = FALSE,
                        n_iter_import = 30000, 
                        sample_every_import = 50, 
                        burnin = 10000)
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

#### measlesoutbreaker analysis toy data: elements of likelihood ####


## Import, no likelihood
data <- outbreaker_data(dates = dt_cases$Date,
                        age_group = dt_cases$age_group,
                        postcode = dt_cases$county,
                        f_dens = matrix(exp(f), nrow = 1),
                        w_dens = NULL, 
                        a_dens = NULL, 
                        genotype = rep("Not attributed", dim(dt_cases)[1]),
                        is_cluster = NULL,
                        population = pop_vect,
                        distance = dist_mat,
                        import = dt_cases$import)

config <- create_config(data = data, 
                        init_tree = "star", 
                        spatial_method = "exponential",
                        init_kappa = 1,
                        gamma = NULL,
                        delta = NULL,
                        n_iter = 70000, 
                        n_iter_import = 30000, 
                        sample_every = 50, 
                        sample_every_import = 50, 
                        burnin = 10000,
                        min_date = -10, 
                        max_kappa = 5,
                        find_import = FALSE,
                        outlier_threshold = 0.05,
                        outlier_relative = FALSE,
                        move_a = FALSE,
                        move_b = FALSE,
                        move_pi = FALSE,
                        move_t_inf = FALSE,
                        move_kappa = FALSE,
                        move_swap_cases = FALSE)
priors <- custom_priors()
likelihoods <- custom_likelihoods(space = f_null, 
                                  timing_infections = f_null,
                                  timing_sampling = f_null,
                                  age = f_null, reporting = f_null )
moves <- custom_moves()
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
saveRDS(out, file = "toy_outbreak_runs/import_no_like.rds")


## Import, time only
data$w_dens = matrix(exp(w), nrow = 1)
data$a_dens = matrix(polymod_prop, nrow = nrow(polymod_prop))
data <- outbreaker_data(data = data)
likelihoods <- custom_likelihoods(space = f_null, 
                                  age = f_null)
config$delta <- 30
config$move_t_inf <- TRUE
config$move_kappa <- TRUE
config$move_swap_cases <- TRUE
config$init_alpha <- NULL
config <- create_config(data = data,
                        config)
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
config$move_a <- FALSE
config$move_b <- FALSE
config$gamma <- NULL
config <- create_config(data = data,
                        config)
likelihoods <- custom_likelihoods(space = f_null)
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
saveRDS(out, file = paste0("toy_outbreak_runs/import_time_genotype_age.rds"))


#### measlesoutbreaker analysis toy data: Different genotype report rate ####

# Set seed
set.seed(1)
# Population per county
pop_county <- dt_distance[!duplicated(county1), pop_county1]
names(pop_county) <- dt_distance[!duplicated(county1), county1]
# Draw the r0 in each state, following an exponential distribution
r0_state <- rexp(n = length(unique(dt_state_county$STNAME)), 
                 rate = 1.5)
names(r0_state) <- unique(dt_state_county$STNAME)
r0_state[r0_state>1] <- 0.95
# Spatial parameters
a <- 0.84
b <- 0.08
# Spatial threshold
gamma <- 100
# Maximum number of cases in the generated dataset
nb_cases <- 1000
dt_cases <- toy_outbreak[["cases"]]
dt_cases <- dt_cases[order(Date), ]

toy_outbreak <- generate_dataset(a = a, b = b, gamma = gamma, 
                                 dt_distance = dt_distance,
                                 polymod_prop = polymod_prop, w = w, 
                                 nb_cases = nb_cases, r0_state = r0_state, 
                                 pop_county = pop_county, 
                                 dt_state_county = dt_state_county,
                                 t_min = as.Date("2010-01-01"), 
                                 t_max = as.Date("2016-01-01"), prop_gen = 1) 

## No import, threshold = 0.05, all genotype reported
data <- outbreaker_data(dates = dt_cases$Date,
                        age_group = dt_cases$age_group,
                        postcode = dt_cases$county,
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
                        n_iter = 70000, 
                        sample_every = 50, 
                        min_date = -10, 
                        max_kappa = 5,
                        find_import = TRUE, 
                        outlier_threshold = 0.05, 
                        outlier_relative = FALSE,
                        n_iter_import = 30000, 
                        sample_every_import = 50, 
                        burnin = 10000)
priors <- custom_priors()
likelihoods <- custom_likelihoods()
moves <- custom_moves()

out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)

saveRDS(out, file = "toy_outbreak_runs/no_import_thresh005_all_gen.rds")


# Set seed
set.seed(1)
# Population per county
pop_county <- dt_distance[!duplicated(county1), pop_county1]
names(pop_county) <- dt_distance[!duplicated(county1), county1]
# Draw the r0 in each state, following an exponential distribution
r0_state <- rexp(n = length(unique(dt_state_county$STNAME)), 
                 rate = 1.5)
names(r0_state) <- unique(dt_state_county$STNAME)
r0_state[r0_state>1] <- 0.95
# Spatial parameters
a <- 0.84
b <- 0.08
# Spatial threshold
gamma <- 100
# Maximum number of cases in the generated dataset
nb_cases <- 1000
dt_cases <- toy_outbreak[["cases"]]
dt_cases <- dt_cases[order(Date), ]

toy_outbreak <- generate_dataset(a = a, b = b, gamma = gamma, 
                                 dt_distance = dt_distance,
                                 polymod_prop = polymod_prop, w = w, 
                                 nb_cases = nb_cases, r0_state = r0_state, 
                                 pop_county = pop_county, 
                                 dt_state_county = dt_state_county,
                                 t_min = as.Date("2010-01-01"), 
                                 t_max = as.Date("2016-01-01"), prop_gen = 0) 

## No import, threshold = 0.05, no genotype reported
data <- outbreaker_data(dates = dt_cases$Date,
                        age_group = dt_cases$age_group,
                        postcode = dt_cases$county,
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
                        n_iter = 70000, 
                        sample_every = 50, 
                        min_date = -10, 
                        max_kappa = 5,
                        find_import = TRUE, 
                        outlier_threshold = 0.05, 
                        outlier_relative = FALSE,
                        n_iter_import = 30000, 
                        sample_every_import = 50, 
                        burnin = 10000)
priors <- custom_priors()
likelihoods <- custom_likelihoods()
moves <- custom_moves()

out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)

saveRDS(out, file = "toy_outbreak_runs/no_import_thresh005_no_gen.rds")
