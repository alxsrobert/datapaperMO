source("R/library_importation.R")
data("fake_outbreak")
source("R/prepare_for_figures.R")

dt_cases <- fake_outbreak[["cases"]]
dist_mat <- fake_outbreak[["distance"]]
pop_vect <- fake_outbreak[["population"]]
age_contact <- fake_outbreak[["age_contact"]]
burnin <- 10000
sample_every <- 50
dt_cases <- dt_cases[order(Date), ]
max_clust <- 500

dt_cases[, size_cluster := table(cluster)[cluster]]

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


#### Reference clusters ####

ref_size_clusters <- dt_cases[!duplicated(cluster) & cluster != ".", size_cluster]
ref_size_clusters <- c(rep(1, dim(dt_cases[cluster == ".",])[1]), ref_size_clusters)
ref_size_clusters <- factor(ref_size_clusters, levels = 1:max_clust)
ref_size_clusters <- table(ref_size_clusters)
size_cluster_ref_barplot <- t(rowsum(x = as.matrix(ref_size_clusters), 
                                     group = groups_barplot))
colnames(size_cluster_ref_barplot) <- names_barplot

size_cluster_ref_singletons <- numeric(length(size_cluster_ref_barplot))
names(size_cluster_ref_singletons) <- names(size_cluster_ref_barplot)
size_cluster_ref_singletons[1] <- size_cluster_ref_barplot[1]

#### threshold 001 ####

out <- readRDS(file = "fake_outbreak_runs/no_import_thresh001.rds")
fig_001 <- prepare_for_figures(out = out, 
                               dt_cases = dt_cases, burnin = burnin, 
                               sample_every = sample_every, 
                               max_clust = max_clust,
                               thresh_barplot = thresh_barplot, diff = diff)

low_size_cluster_barplot_001 <- fig_001[["low_size_cluster_barplot"]]
med_size_cluster_barplot_001 <- fig_001[["med_size_cluster_barplot"]]
up_size_cluster_barplot_001 <- fig_001[["up_size_cluster_barplot"]]
med_prop_singletons_001 <- fig_001[["med_prop_singletons"]]
low_001 <- fig_001[["low"]]
med_001 <- fig_001[["med"]]
up_001 <- fig_001[["up"]]
dt_heatmap_001 <- fig_001[["dt_heatmap"]]
dt_heatmap_001[, type := "001"]

#### threshold 005 ####

out <- readRDS(file = "fake_outbreak_runs/no_import_thresh005.rds")

fig_005 <- prepare_for_figures(out = out, 
                               dt_cases = dt_cases, burnin = burnin, 
                               sample_every = sample_every, 
                               max_clust = max_clust,
                               thresh_barplot = thresh_barplot, diff = diff)

low_size_cluster_barplot_005 <- fig_005[["low_size_cluster_barplot"]]
med_size_cluster_barplot_005 <- fig_005[["med_size_cluster_barplot"]]
up_size_cluster_barplot_005 <- fig_005[["up_size_cluster_barplot"]]
med_prop_singletons_005 <- fig_005[["med_prop_singletons"]]
low_005 <- fig_005[["low"]]
med_005 <- fig_005[["med"]]
up_005 <- fig_005[["up"]]
dt_heatmap_005 <- fig_005[["dt_heatmap"]]
dt_heatmap_005[, type := "005"]

#### threshold 09 ####

out <- readRDS(file = "fake_outbreak_runs/no_import_thresh09.rds")

fig_09 <- prepare_for_figures(out = out, 
                              dt_cases = dt_cases, burnin = burnin, 
                              sample_every = sample_every, 
                              max_clust = max_clust,
                              thresh_barplot = thresh_barplot, diff = diff)

low_size_cluster_barplot_09 <- fig_09[["low_size_cluster_barplot"]]
med_size_cluster_barplot_09 <- fig_09[["med_size_cluster_barplot"]]
up_size_cluster_barplot_09 <- fig_09[["up_size_cluster_barplot"]]
med_prop_singletons_09 <- fig_09[["med_prop_singletons"]]
low_09 <- fig_09[["low"]]
med_09 <- fig_09[["med"]]
up_09 <- fig_09[["up"]]
dt_heatmap_09 <- fig_09[["dt_heatmap"]]
dt_heatmap_09[, type := "09"]

#### threshold 095 ####

out <- readRDS(file = "fake_outbreak_runs/no_import_thresh095.rds")

fig_095 <- prepare_for_figures(out = out, 
                               dt_cases = dt_cases, burnin = burnin, 
                               sample_every = sample_every, 
                               max_clust = max_clust,
                               thresh_barplot = thresh_barplot, diff = diff)

low_size_cluster_barplot_095 <- fig_095[["low_size_cluster_barplot"]]
med_size_cluster_barplot_095 <- fig_095[["med_size_cluster_barplot"]]
up_size_cluster_barplot_095 <- fig_095[["up_size_cluster_barplot"]]
med_prop_singletons_095 <- fig_095[["med_prop_singletons"]]
low_095 <- fig_095[["low"]]
med_095 <- fig_095[["med"]]
up_095 <- fig_095[["up"]]
dt_heatmap_095 <- fig_095[["dt_heatmap"]]
dt_heatmap_095[, type := "095"]

#### with imports ####

out <- readRDS(file = "fake_outbreak_runs/with_import.rds")

fig_import <- prepare_for_figures(out = out, 
                                  dt_cases = dt_cases, burnin = burnin, 
                                  sample_every = sample_every, 
                                  max_clust = max_clust,
                                  thresh_barplot = thresh_barplot, diff = diff)

low_size_cluster_barplot_import <- fig_import[["low_size_cluster_barplot"]]
med_size_cluster_barplot_import <- fig_import[["med_size_cluster_barplot"]]
up_size_cluster_barplot_import <- fig_import[["up_size_cluster_barplot"]]
med_prop_singletons_import <- fig_import[["med_prop_singletons"]]
low_import <- fig_import[["low"]]
med_import <- fig_import[["med"]]
up_import <- fig_import[["up"]]
dt_heatmap_import <- fig_import[["dt_heatmap"]]
dt_heatmap_import[, type := "import"]

#### with imports 005 ####

out <- readRDS(file = "fake_outbreak_runs/with_import_005.rds")

fig_005_wi <- prepare_for_figures(out = out, 
                                  dt_cases = dt_cases, burnin = burnin, 
                                  sample_every = sample_every, 
                                  max_clust = max_clust,
                                  thresh_barplot = thresh_barplot, diff = diff)

low_size_cluster_barplot_005_wi <- fig_005_wi[["low_size_cluster_barplot"]]
med_size_cluster_barplot_005_wi <- fig_005_wi[["med_size_cluster_barplot"]]
up_size_cluster_barplot_005_wi <- fig_005_wi[["up_size_cluster_barplot"]]
med_prop_singletons_005_wi <- fig_005_wi[["med_prop_singletons"]]
low_005_wi <- fig_005_wi[["low"]]
med_005_wi <- fig_005_wi[["med"]]
up_005_wi <- fig_005_wi[["up"]]
dt_heatmap_005_wi <- fig_005_wi[["dt_heatmap"]]
dt_heatmap_005_wi[, type := "005_wi"]


#### with imports 095 ####

out <- readRDS(file = "fake_outbreak_runs/with_import_095.rds")

fig_095_wi <- prepare_for_figures(out = out, 
                                  dt_cases = dt_cases, burnin = burnin, 
                                  sample_every = sample_every, 
                                  max_clust = max_clust,
                                  thresh_barplot = thresh_barplot, diff = diff)

low_size_cluster_barplot_095_wi <- fig_095_wi[["low_size_cluster_barplot"]]
med_size_cluster_barplot_095_wi <- fig_095_wi[["med_size_cluster_barplot"]]
up_size_cluster_barplot_095_wi <- fig_095_wi[["up_size_cluster_barplot"]]
med_prop_singletons_095_wi <- fig_095_wi[["med_prop_singletons"]]
low_095_wi <- fig_095_wi[["low"]]
med_095_wi <- fig_095_wi[["med"]]
up_095_wi <- fig_095_wi[["up"]]
dt_heatmap_095_wi <- fig_095_wi[["dt_heatmap"]]
dt_heatmap_095_wi[, type := "095_wi"]