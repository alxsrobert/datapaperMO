## Generate all the figures (main and supplement) from analysis of toy_outbreak

# Load all runs
source("R/load_all_analysis.R")

#### Generate figure ####

## Call function to generate figures 3
ref_breaks <- c(1, 10, 50, 100, 500)
generate_figure_3(dt_cases, ref_breaks)

## Call function to generate figures 4 and 5
# NOTE: In the simulations there is only 1 "reference" import status distribution, 
# compared with the "epi" and "ideal" distributions in the US dataset. 
fig_hist_list <- list(fig_import, fig_ref)
list_fig_heatmap <- list(fig_import)

generate_figure_4_5(fig_hist_list, list_fig_heatmap)
par(fig = c(0.5, 1, 0.5, 1))
legend("topright", fill = c(grey.colors(length(fig_hist_list)), NA), cex = 1.2,
       border = NA,  pch = c(NA, NA, "_"), pt.cex = 3,
       legend = c("Epi import", "Epidemiological cluster",
                  "Correct inference"), bty = "n")

fig_hist_list <- list(fig_001, fig_005, fig_095, fig_005_wi, fig_ref)
list_fig_heatmap <- list(fig_005, fig_005_wi)

generate_figure_4_5(fig_hist_list, list_fig_heatmap)
par(fig = c(0.5, 1, 0.5, 1))
legend("right", fill = c(grey.colors(length(fig_hist_list)), NA), 
       cex = 1.2, border = NA, pch = c(NA, NA, NA, NA, NA, "_"), pt.cex = 3,
       legend = c("Absolute threshold 1%", "Absolute Threshold 5%", 
                  "Relative threshold 95%","Absolute threshold 5% + imports",
                  "Epidemiological cluster", "Correct inference"), bty = "n")
## Subsequent cases / number of imports in each state at every iteration
list_factor_import <- list(fig_005[["factor_import"]], 
                           fig_005_wi[["factor_import"]])
## Categories in figure 6
ref_breaks <- c(0, 1, 3, 5, 10, 50)

generate_figure_6(list_factor_import, ref_breaks)

#### Plot supplement ####

## Impact of the proportion of genotyped cases
fig_hist_list <- list(fig_no_gen_sim, fig_005_sim, fig_all_gen_sim, fig_sim)
list_fig_heatmap <- list(fig_no_gen_sim, fig_all_gen_sim)

generate_figure_4_5(fig_hist_list, list_fig_heatmap)
par(fig = c(0.5, 1, 0.5, 1))
legend("right", fill = c(grey.colors(length(fig_hist_list)), NA), 
       cex = 1.2, border = NA, pch = c(NA, NA, NA, NA, NA, "_"), pt.cex = 3,
       legend = c("No case genotyped", "60% cases genotyped","All cases genotyped",
                  "Epidemiological cluster", "Correct inference"), bty = "n")

##Section 1: description data
supp_desc_data(dt_cases = dt_cases)


## Precision and sensitivity clustering of cases
list_fig <- list(fig_import, fig_005, fig_005_wi, 
                 fig_09, fig_095)
titles <- c("With import", "Absolute 0.05", 
            "Absolute 0.05 with import", "Relative 90%", "Relative 95%")
supp_fig_eval_cluster(list_fig, titles = titles)

## Posterior traces
list_out <- list(out_import, out_005, out_005_wi, 
                 out_09, out_095)
supp_post(list_out, burnin, sample_every)

## cluster stratified by state
supp_fig_stratified_state(out_import, dt_cases, burnin, sample_every, max_clust)

supp_fig_stratified_state(out_005, dt_cases, burnin, sample_every, max_clust)

supp_fig_stratified_state(out_095, dt_cases, burnin, sample_every, max_clust)

supp_fig_stratified_state(out_005_wi, dt_cases, burnin, sample_every, max_clust)

## Parameter estimates
list_out <- list(out_import, out_005, out_09,
                 out_095, out_005_wi, out_095_wi)

supp_fig_param_estimate(list_out, burnin, sample_every)
mgp.axis(side = 1, at = 1:length(list_out), 
         labels = c("With import\n", "Absolute 0.05\n", "Relative 90%\n", 
                    "Relative 95%\n","Absolute 0.05 \nw/ imports",
                    "Relative 95% \nw/ imports"), 
         mgp = c(3.5,2.5,0))

## Distance transmission
list_out <- list(out_import, out_005, out_005_wi)
supp_fig_distance_transmission(list_out = list_out, burnin = burnin, 
                               sample_every = sample_every, dt_cases = dt_cases)
legend("right", fill = grey.colors(length(list_out)), cex = 1.4, border = NA, 
       legend = c("With import", "Absolute threshold 5%", 
                  "Absolute threshold 5% + imports"), bty = "n")

## Sensitivity components of likelihood 
fig_hist_list <- list(fig_time, fig_time_gen, fig_time_age, fig_time_spa, fig_ref)
list_fig_heatmap <- list(fig_time_gen, fig_time_spa)

generate_figure_4_5(fig_hist_list, list_fig_heatmap)
par(fig = c(0.5, 1, 0.5, 1))
legend("right", fill = c(grey.colors(length(fig_hist_list)), NA), cex = 1.2,
       border = NA,  pch = c(rep(NA, length(fig_hist_list)), "_"), pt.cex = 3,
       legend = c("Time only", "Time and genotype", "Time, genotype and age",
                  "Time, genotype and space", "Epidemiological cluster",
                  "Correct inference"),
       bty = "n")

## Number of secondary transmission:
# Overall
list_out <- list(out_import, out_005, out_005_wi)
supp_fig_sec_overall(list_out, burnin, sample_every)
legend("right", fill = grey.colors(length(list_out)), cex = 1.4, border = NA, 
       legend = c("With import", "Absolute threshold 5%", 
                  "Absolute threshold 5% + imports"), bty = "n")
# Stratified by state
list_out <- list(out_import, out_005, out_005_wi)
supp_fig_sec_map(list_out, dt_cases)
