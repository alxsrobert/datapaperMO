
library(data.table)
library(measlesoutbreaker)
library(geosphere)
library(xlsx)
library(dplyr)
library(Rcpp)


data("fake_outbreak")

dt_cases <- fake_outbreak[["cases"]]
dist_mat <- fake_outbreak[["distance"]]
pop_vect <- fake_outbreak[["population"]]
age_contact <- fake_outbreak[["age_contact"]]

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
                        outlier_threshold = 0.05, 
                        outlier_relative = FALSE,
                        n_iter_import = 30000, 
                        sample_every_import = 50, 
                        burnin = 10000)


priors <- custom_priors()
likelihoods <- custom_likelihoods()
moves <- custom_moves()
start <- Sys.time()
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
end <- Sys.time()
(end-start) %>% print
# Time difference of  hours

saveRDS(out, file = "no_import_thresh005.rds")



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
                        import = dt_cases$import)

config <- create_config(data = data, 
                        init_tree = "star", 
                        spatial_method = "exponential",
                        gamma = 100,
                        delta = 30,
                        n_iter = 70000, 
                        sample_every = 50, 
                        min_date = -10, 
                        max_kappa = 5,
                        find_import = FALSE, 
                        outlier_threshold = 0.05, 
                        outlier_relative = FALSE,
                        n_iter_import = 30000, 
                        sample_every_import = 50, 
                        burnin = 10000)


priors <- custom_priors()
likelihoods <- custom_likelihoods()
moves <- custom_moves()
start <- Sys.time()
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
end <- Sys.time()
(end-start) %>% print
# Time difference of 1.66 hours
saveRDS(out, file = "with_import.rds")

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
start <- Sys.time()
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
end <- Sys.time()
(end-start) %>% print
# Time difference of 2.43 hours

saveRDS(out, file = "no_import_thresh001.rds")

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
                        outlier_threshold = 0.95, 
                        outlier_relative = TRUE,
                        n_iter_import = 30000, 
                        sample_every_import = 50, 
                        burnin = 10000)


priors <- custom_priors()
likelihoods <- custom_likelihoods()
moves <- custom_moves()
start <- Sys.time()
out <- outbreaker(data = data,config = config,priors = priors,
                  likelihoods = likelihoods,moves = moves)
end <- Sys.time()
(end-start) %>% print
# Time difference of 2.671993 hours

saveRDS(out, file = "no_import_thresh095.rds")
