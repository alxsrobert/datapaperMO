source("R/load_analysis_likelihoods.R")
source("R/load_analysis_for_figure.R")

#### Generate figure ####

## Call function to generate figures 3
ref_breaks <- c(1, 10, 50, 100, 500)
categ <- c("1-10", "10-50", "50-100", "100-500")
generate_figure_3(dt_cases, ref_breaks, categ)

## Call function to generate figures 4 and 5
fig_hist_list <- list(fig_import, fig_ref)
list_fig_heatmap <- list(fig_import)

generate_figure_4_5(fig_hist_list, list_fig_heatmap)

fig_hist_list <- list(fig_001, fig_005, fig_095, fig_005_wi, fig_ref)
list_fig_heatmap <- list(fig_005, fig_005_wi)

generate_figure_4_5(fig_hist_list, list_fig_heatmap)

## Subsequent cases / number of imports in each state at every iteration
list_factor_import <- list(fig_005[["factor_import"]], 
                           fig_005_wi[["factor_import"]])
## Categories in figure 5
ref_breaks <- c(0, 1, 3, 5, 10, 50)
categ <- c("0-1",
           "1-3",
           "3-5",
           "5-10",
           "10+")

generate_figure_6(list_factor_import, ref_breaks, categ)

#### Plot supplement ####

list_out <- list(out_import, out_005, out_09,
                 out_095, out_005_wi, out_095_wi)

supp_fig_param_estimate(list_out, burnin, sample_every)

supp_fig_stratified_state(out_import, dt_cases, burnin, sample_every, max_clust)

supp_fig_stratified_state(out_005, dt_cases, burnin, sample_every, max_clust)

supp_fig_stratified_state(out_095, dt_cases, burnin, sample_every, max_clust)

supp_fig_stratified_state(out_005_wi, dt_cases, burnin, sample_every, max_clust)

list_out <- list(out_import, out_005, out_005_wi)
supp_fig_sec_overall(list_out, burnin, sample_every)

list_out <- list(out_ideal, out_005, out_005_wi)
supp_fig_sec_map(list_out, dt_cases)

list_out <- list(out_ideal, out_005, out_005_wi)
supp_fig_distance_transmission(list_out = list_out, burnin = burnin, 
                               sample_every = sample_every, dt_cases = dt_cases,
                               dt_distance = dt_distance)

list_fig <- list(fig_import, fig_005, fig_005_wi, 
                 fig_09, fig_095)
supp_fig_eval_cluster(list_fig, titles = NULL)

list_out <- list(out_import, out_005, out_005_wi, 
                 out_09, out_095)
supp_post(list_out, burnin, sample_every)


supp_desc_data(dt_cases = dt_cases)

## Sensitivity

fig_hist_list <- list(fig_no_like, fig_time, fig_time_gen, fig_time_age, fig_time_spa, fig_ref)
list_fig_heatmap <- list(fig_time_gen, fig_time_spa)

generate_figure_4_5(fig_hist_list, list_fig_heatmap)

fig_hist_list <- list(fig_no_like, fig_time, fig_time_gen, fig_time_age, fig_time_spa)
generate_figure_4_5(fig_hist_list, list_fig_heatmap, ref = F)
