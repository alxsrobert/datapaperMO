#### Load libraries and parameters ####

source("R/library_importation.R")
data("toy_outbreak")
source("R/function_prepare_for_figures.R")
source("R/function_generate_figures_main.R")
source("R/function_supplement_figures.R")
## Load elements of the list toy_outbreak
dt_cases <- toy_outbreak[["cases"]]
dt_cases <- dt_cases[order(Date), ]
dt_cases[, size_cluster := table(cluster)[cluster]]
dist_mat <- toy_outbreak[["distance"]]
pop_vect <- toy_outbreak[["population"]]
age_contact <- toy_outbreak[["age_contact"]]

burnin <- 10000
sample_every <- 50
max_clust <- 500

## Break points and labels for histogram
thresh_barplot <- c(1,2,3,5,10,15,20,30,40,50)
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

# Gather the reference values in the list fig_sim
fig_sim <- list(med_size_cluster_barplot = size_cluster_ref_barplot,
                med_prop_singletons = size_cluster_ref_barplot,
                med = med_ref)

#### Load measlesoutbreaker runs ####

## Inference with absolute threshold (k = 0.05) when no cases was genotyped
out_no_gen_sim <- readRDS(file = "toy_outbreak_runs/no_import_thresh005_no_gen.rds")
fig_no_gen_sim <- prepare_for_figures(out = out_no_gen_sim,
                                      dt_cases = dt_cases, burnin = burnin,
                                      sample_every = sample_every,
                                      max_clust = max_clust,
                                      thresh_barplot = thresh_barplot, diff = diff)

## Inference with absolute threshold (k = 0.05) when all cases were genotyped
out_all_gen_sim <- readRDS(file = "toy_outbreak_runs/no_import_thresh005_all_gen.rds")
fig_all_gen_sim <- prepare_for_figures(out = out_all_gen_sim, dt_cases = dt_cases, burnin = burnin,
                                       sample_every = sample_every, max_clust = max_clust,
                                       thresh_barplot = thresh_barplot, diff = diff)

## Inference with absolute threshold (k = 0.05) and 40% of genotyped cases (same as data)
out_005_sim <- readRDS(file = "toy_outbreak_runs/no_import_thresh005.rds")
fig_005_sim <- prepare_for_figures(out = out_005_sim, dt_cases = dt_cases, burnin = burnin,
                                   sample_every = sample_every, max_clust = max_clust,
                                   thresh_barplot = thresh_barplot, diff = diff)

