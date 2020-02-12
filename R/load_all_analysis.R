## Load o2geosocial runs comparing the impact of the import threshold lambda 
#### Load libraries and parameters ####

source("R/library_importation.R")
data("toy_outbreak_long")
source("R/function_prepare_for_figures.R")
source("R/function_generate_figures_main.R")
source("R/function_supplement_figures.R")
## Load elements of the list toy_outbreak_long
dt_cases <- toy_outbreak_long[["cases"]]
dt_cases <- dt_cases[order(Date), ]
dt_cases[, size_cluster := table(cluster)[as.character(cluster)]]
dt_regions <- toy_outbreak_long[["dt_regions"]]
all_dist <- distGeo(matrix(c(rep(dt_regions$long, nrow(dt_regions)), 
                             rep(dt_regions$lat, nrow(dt_regions))), ncol = 2), 
                    matrix(c(rep(dt_regions$long, each = nrow(dt_regions)), 
                             rep(dt_regions$lat, each = nrow(dt_regions))), ncol = 2))

dist_mat <- matrix(all_dist/1000, nrow = nrow(dt_regions))
pop_vect <- dt_regions$population
age_contact <- toy_outbreak_long[["age_contact"]]

burnin <- 10000
sample_every <- 50
max_clust <- 500

## Break points and labels for histogram
thresh_barplot <- c(1,2,3,5,10, 20,30,50, 100, 200)
names_barplot <- sapply(1:length(thresh_barplot), function(X){
  if((X + 1) > length(thresh_barplot))
    return(paste0(as.character(thresh_barplot[X]), "+"))
  if(thresh_barplot[X+1] == (thresh_barplot[X] + 1))
    return(as.character(thresh_barplot[X]))
  else 
    return(paste0(as.character(thresh_barplot[X]), "-",
                  as.character(thresh_barplot[X+1] - 1)))
})
groups_barplot <- sapply(1:max_clust, function(X){
  return(max(which(thresh_barplot <= X)))
})
diff <- 0.25

## Load data with all abbreviations and state names
data(state, package = "datasets")
dt_map_cases <- data.table(abb = c(state.abb, "COL"))
dt_map_cases[, region := tolower(c(state.name, "district of columbia"))]
setkey(dt_map_cases, region)
dt_map_cases[, cases := 0]
dt_cases$State <- factor(x = dt_map_cases[tolower(dt_cases$State), abb], 
                         levels = dt_map_cases$abb)


#### Reference clusters ####

## Size of each cluster in data
ref_size_clusters <- dt_cases[!duplicated(cluster) & cluster != ".", size_cluster]
ref_size_clusters <- factor(ref_size_clusters, levels = 1:max_clust)
ref_size_clusters <- table(ref_size_clusters)

## Number of cluster in each group from group_barplot
size_cluster_ref_barplot <- t(rowsum(x = as.matrix(ref_size_clusters), 
                                     group = groups_barplot)) %>% as.numeric
names(size_cluster_ref_barplot) <- names_barplot

## Number of singletons in data
size_cluster_ref_singletons <- numeric(length(size_cluster_ref_barplot))
names(size_cluster_ref_singletons) <- names(size_cluster_ref_barplot)
size_cluster_ref_singletons[1] <- size_cluster_ref_barplot[1]

med_ref <- rep(size_cluster_ref_barplot %>% sum, 2)
names(med_ref) <- c("Imports", "Import status \ncorrectly inferred")

fig_sim <- fig_ref <- list(med_size_cluster_barplot = size_cluster_ref_barplot,
                           med_prop_singletons = size_cluster_ref_barplot,
                           med = med_ref)

#### Load all analysis ####

## threshold 001
out_001 <- readRDS(file = "toy_outbreak_runs/no_import_thresh001.rds")
fig_001 <- readRDS(file = "fig_runs/fig_001.rds")
## threshold 005
out_005 <- readRDS(file = "toy_outbreak_runs/no_import_thresh005.rds")
fig_005 <- readRDS(file = "fig_runs/fig_005.rds")
## threshold 09
out_09 <- readRDS(file = "toy_outbreak_runs/no_import_thresh09.rds")
fig_09 <- readRDS(file = "fig_runs/fig_09.rds")
## threshold 095
out_095 <- readRDS(file = "toy_outbreak_runs/no_import_thresh095.rds")
fig_095 <- readRDS(file = "fig_runs/fig_095.rds")
## with imports
out_import <- readRDS(file = "toy_outbreak_runs/with_import.rds")
fig_import <- readRDS(file = "fig_runs/fig_import.rds")
## with imports 005
out_005_wi <- readRDS(file = "toy_outbreak_runs/with_import_005.rds")
fig_005_wi <- readRDS(file = "fig_runs/fig_005_wi.rds")
## with imports 095
out_095_wi <- readRDS(file = "toy_outbreak_runs/with_import_095.rds")
fig_095_wi <- readRDS(file = "fig_runs/fig_095_wi.rds")
## time only
out_time <- readRDS(file = "toy_outbreak_runs/import_time_only.rds")
fig_time <- readRDS(file = "fig_runs/fig_time.rds")
## time and genotype
out_time_gen <- readRDS(file = "toy_outbreak_runs/import_time_genotype.rds")
fig_time_gen <- readRDS(file = "fig_runs/fig_time_gen.rds")
## time, genotype and age
out_time_age <- readRDS(file = "toy_outbreak_runs/import_time_genotype_age.rds")
fig_time_age <- readRDS(file = "fig_runs/fig_time_age.rds")
## time, genotype and sex
out_time_spa <- readRDS(file = "toy_outbreak_runs/import_time_genotype_space.rds")
fig_time_spa <- readRDS(file = "fig_runs/fig_time_spa.rds")

## Inference with absolute threshold (k = 0.05) when no cases was genotyped
out_no_gen_sim <- readRDS(file = "toy_outbreak_runs/no_import_thresh005_no_gen.rds")
fig_no_gen_sim <- readRDS(file = "fig_runs/fig_no_gen_sim.rds")

## Inference with absolute threshold (k = 0.05) when all cases were genotyped
out_all_gen_sim <- readRDS(file = "toy_outbreak_runs/no_import_thresh005_all_gen.rds")
fig_all_gen_sim <- readRDS(file = "fig_runs/fig_all_gen_sim.rds")

## Inference with absolute threshold (k = 0.05) and 40% of genotyped cases (same as data)
out_005_sim <- readRDS(file = "toy_outbreak_runs/no_import_thresh005.rds")
fig_005_sim <- readRDS(file = "fig_runs/fig_005_sim.rds")

