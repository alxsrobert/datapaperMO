#### Load libraries and parameters ####

source("R/library_importation.R")
data("toy_outbreak")
source("R/prepare_for_figures.R")
source("R/function_generate_figures_main.R")
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
                                     group = groups_barplot))
colnames(size_cluster_ref_barplot) <- names_barplot

## Number of singletons in data
size_cluster_ref_singletons <- numeric(length(size_cluster_ref_barplot))
names(size_cluster_ref_singletons) <- names(size_cluster_ref_barplot)
size_cluster_ref_singletons[1] <- size_cluster_ref_barplot[1]

#### Load measlesoutbreaker runs ####

## threshold 001
out <- readRDS(file = "toy_outbreak_runs/no_import_thresh001.rds")
fig_001 <- prepare_for_figures(out = out, 
                               dt_cases = dt_cases, burnin = burnin, 
                               sample_every = sample_every, 
                               max_clust = max_clust,
                               thresh_barplot = thresh_barplot, diff = diff)

## threshold 005
out <- readRDS(file = "toy_outbreak_runs/no_import_thresh005.rds")
fig_005 <- prepare_for_figures(out = out, dt_cases = dt_cases, burnin = burnin, 
                               sample_every = sample_every, max_clust = max_clust,
                               thresh_barplot = thresh_barplot, diff = diff)

## threshold 09
out <- readRDS(file = "toy_outbreak_runs/no_import_thresh09.rds")
fig_09 <- prepare_for_figures(out = out, 
                              dt_cases = dt_cases, burnin = burnin, 
                              sample_every = sample_every, 
                              max_clust = max_clust,
                              thresh_barplot = thresh_barplot, diff = diff)

## threshold 095
out <- readRDS(file = "toy_outbreak_runs/no_import_thresh095.rds")
fig_095 <- prepare_for_figures(out = out, 
                               dt_cases = dt_cases, burnin = burnin, 
                               sample_every = sample_every, 
                               max_clust = max_clust,
                               thresh_barplot = thresh_barplot, diff = diff)

## with imports
out <- readRDS(file = "toy_outbreak_runs/with_import.rds")
fig_import <- prepare_for_figures(out = out, 
                                  dt_cases = dt_cases, burnin = burnin, 
                                  sample_every = sample_every, 
                                  max_clust = max_clust,
                                  thresh_barplot = thresh_barplot, diff = diff)

## with imports 005
out <- readRDS(file = "toy_outbreak_runs/with_import_005.rds")
fig_005_wi <- prepare_for_figures(out = out, 
                                  dt_cases = dt_cases, burnin = burnin, 
                                  sample_every = sample_every, 
                                  max_clust = max_clust,
                                  thresh_barplot = thresh_barplot, diff = diff)

## with imports 095
out <- readRDS(file = "toy_outbreak_runs/with_import_095.rds")
fig_095_wi <- prepare_for_figures(out = out, 
                                  dt_cases = dt_cases, burnin = burnin, 
                                  sample_every = sample_every, 
                                  max_clust = max_clust,
                                  thresh_barplot = thresh_barplot, diff = diff)

#### Generate figure ####

## Add id to heatmaps
fig_001[["dt_heatmap"]][, type := "001"]
fig_005[["dt_heatmap"]][, type := "005"]
fig_09[["dt_heatmap"]][, type := "09"]
fig_095[["dt_heatmap"]][, type := "095"]
fig_import[["dt_heatmap"]][, type := "import"]
fig_005_wi[["dt_heatmap"]][, type := "005_wi"]

## List of all heatmaps in data
list_heatmap <- list(#fig_001[["dt_heatmap"]],
                     fig_005[["dt_heatmap"]],
                     # fig_09[["dt_heatmap"]],
                     # fig_095[["dt_heatmap"]], 
                     # fig_import[["dt_heatmap"]],
                     fig_005_wi[["dt_heatmap"]]
                     )
## Median values of cluster size distributions
med_size_cluster <- rbind(fig_001[["med_size_cluster_barplot"]],
                          fig_005[["med_size_cluster_barplot"]],
                          fig_09[["med_size_cluster_barplot"]],
                          fig_095[["med_size_cluster_barplot"]],
                          fig_import[["med_size_cluster_barplot"]],
                          fig_005_wi[["med_size_cluster_barplot"]],
                          size_cluster_ref_barplot)
## 97.5% CI
up_size_cluster <- rbind(fig_001[["up_size_cluster_barplot"]],
                         fig_005[["up_size_cluster_barplot"]],
                         fig_09[["up_size_cluster_barplot"]],
                         fig_095[["up_size_cluster_barplot"]],
                         fig_import[["up_size_cluster_barplot"]],
                         fig_005_wi[["up_size_cluster_barplot"]])
## 2.5% CI
low_size_cluster <- rbind(fig_001[["low_size_cluster_barplot"]],
                          fig_005[["low_size_cluster_barplot"]],
                          fig_09[["low_size_cluster_barplot"]],
                          fig_095[["low_size_cluster_barplot"]],
                          fig_import[["low_size_cluster_barplot"]],
                          fig_005_wi[["low_size_cluster_barplot"]])
## Number of singletons
singletons <- rbind(fig_001[["med_prop_singletons"]],
                    fig_005[["med_prop_singletons"]],
                    fig_09[["med_prop_singletons"]],
                    fig_095[["med_prop_singletons"]],
                    fig_import[["med_prop_singletons"]],
                    fig_005_wi[["med_prop_singletons"]],
                    size_cluster_ref_singletons)
## Median number of imports, and number of imports correctly inferred
med_imports<- rbind(fig_001[["med"]], 
                    fig_005[["med"]], 
                    fig_09[["med"]], 
                    fig_095[["med"]], 
                    fig_import[["med"]], 
                    fig_005_wi[["med"]],
                    size_cluster_ref_barplot %>% sum)
## 97.5% CI
up_imports<- rbind(fig_001[["up"]],
                   fig_005[["up"]],
                   fig_09[["up"]],
                   fig_095[["up"]],
                   fig_import[["up"]],
                   fig_005_wi[["up"]])
## 2.5% CI
low_imports<- rbind(fig_001[["low"]],
                    fig_005[["low"]],
                    fig_09[["low"]],
                    fig_095[["low"]],
                    fig_import[["low"]],
                    fig_005_wi[["low"]])
## Call function to generate figures 3 and 4
generate_figure_3_4(med_size_cluster,
                    up_size_cluster,
                    low_size_cluster,
                    singletons,
                    med_imports,
                    up_imports,
                    low_imports,
                    list_heatmap)

## Subsequent cases / number of imports in each state at every iteration
list_factor_import <- list(fig_005[["factor_import"]], 
                           fig_09[["factor_import"]],
                           fig_005_wi[["factor_import"]])
## Categories in figure 5
ref_breaks <- c(0, 1, 3, 5, 10, 50)
categ <- c("0-1",
           "1-3",
           "3-5",
           "5-10",
           "10+")

generate_figure_5(list_factor_import, ref_breaks, categ)
