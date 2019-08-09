source("R/library_importation.R")
source("R/generate_dataset.R")
source("R/load_data_distributions.R")

generate <- TRUE

if(generate){
  set.seed(1)
  pop_county <- dt_distance[!duplicated(county1), pop_county1]
  names(pop_county) <- dt_distance[!duplicated(county1), county1]
  
  r0_state <- rexp(n = length(unique(dt_state_county$STNAME)), 
                   rate = 1.5)
  names(r0_state) <- unique(dt_state_county$STNAME)
  r0_state[r0_state>1] <- 0.95
  
  a <- 0.84
  b <- 0.08
  gamma <- 100
  nb_cases <- 500
  
  toy_outbreak <- generate_dataset(a = a, b = b, gamma = gamma, 
                                   dt_distance = dt_distance,
                                   dt_population = dt_population,
                                   polymod_prop = polymod_prop, w = w, 
                                   nb_cases = nb_cases, r0_state = r0_state, 
                                   pop_county = pop_county, 
                                   dt_state_county = dt_state_county,
                                   t_min = as.Date("2010-01-01"), 
                                   t_max = as.Date("2016-01-01")) 
} else
  data("toy_outbreak")


#### Analysis toy data ####


dt_cases <- toy_outbreak[["cases"]]
dist_mat <- toy_outbreak[["distance"]]
pop_vect <- toy_outbreak[["population"]]
age_contact <- toy_outbreak[["age_contact"]]

dt_cases <- dt_cases[order(Date), ]


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

saveRDS(out, file = "fake_outbreak_runs/no_import_thresh001.rds")

config$outlier_threshold <- 0.05

config <- create_config(config, data = data)

out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)

saveRDS(out, file = "fake_outbreak_runs/no_import_thresh005.rds")

config$outlier_threshold <- 0.9
config$outlier_relative <- TRUE

config <- create_config(data = data, 
                        config = config)


out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)

saveRDS(out, file = "fake_outbreak_runs/no_import_thresh09.rds")

config$outlier_threshold <- 0.95

config <- create_config(data = data, 
                        config = config)


out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)

saveRDS(out, file = "fake_outbreak_runs/no_import_thresh095.rds")

data$import <- dt_cases$import
config$find_import <- FALSE

data <- outbreaker_data(data = data)
config <- create_config(data = data, 
                        config = config)


out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
saveRDS(out, file = "fake_outbreak_runs/with_import.rds")

config$find_import <- TRUE
config$outlier_threshold <- 0.05
config$outlier_relative <- FALSE

config <- create_config(data = data, 
                        config = config)


out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
saveRDS(out, file = "fake_outbreak_runs/with_import_005.rds")

data <- outbreaker_data(data = data,
                        import = dt_cases$import)

config$outlier_threshold <- 0.95
config$outlier_relative <- TRUE

config <- create_config(data = data, 
                        config = config,
                        outlier_threshold = 0.95, 
                        outlier_relative = TRUE)


out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
saveRDS(out, file = "fake_outbreak_runs/with_import_095.rds")
