# datapaperMO
analysis of toy_outbreak dataset from the measlesoutbreaker package

install package measlesoutbreaker from "alxsrobert/measlesoutbreaker"
```R
install.packages("devtools")
library(devtools)
install.github("alxsrobert/measlesoutbreaker")
library(measlesoutbreaker)
```

In data: 
- State_initials.csv	: Data table with state names and acronyms.
- pop_center.csv: Population centroid and population of every county (Source: https://www.census.gov/geographies/reference-files/2010/geo/2010-centers-population.html (per state : UNITED STATES)).
- social_contact_uk.csv : Contact Matrix of contacts between age groups (Source: Social Contacts and Mixing Patterns Relevant to the Spread of Infectious Diseases)

In R:
- library_importation.R: import libraries for analysis.
- function_generate_dataset.R: Functions to generate toy_outbreak.
- function_prepare_for_figures.R: Function to generate summary statistics on run.
- function_generate_figures_main.R: function to generate figures similar to figure 3, 4, 5 and 6.
- function_supplement_figures.R: function to generate figures similar to figures in the supplement.
- load_analysis_threshold.R: Load measlesoutbreaker runs generated with different import thresholds.
- load_analysis_likelihoods.R: Load measlesoutbreaker runs generated with different likelihoods.
- load_analysis_genotype.R: Load measlesoutbreaker runs generated with different proportion of genotyped cases.
- load_data_distributions.R: Load prior distributions.
- analysis_generated_data.R: Script to run measlesoutbreaker on toy_outbreak.
- generate_figure.R: Generate all figures (main + supplement).

To generate all figures, run the script generate_figure.R.
Script to generate figures 3, 4, 5 and 6:

```R
source("R/load_analysis_threshold.R")
source("R/load_analysis_likelihoods.R")
source("R/load_analysis_genotype.R")
## Call function to generate figures 3
ref_breaks <- c(1, 10, 50, 100, 500)
generate_figure_3(dt_cases, ref_breaks)

## Call function to generate figures 4 and 5
fig_hist_list <- list(fig_import, fig_ref)
list_fig_heatmap <- list(fig_import)

generate_figure_4_5(fig_hist_list, list_fig_heatmap)

fig_hist_list <- list(fig_001, fig_005, fig_095, fig_005_wi, fig_ref)
list_fig_heatmap <- list(fig_005, fig_005_wi)

generate_figure_4_5(fig_hist_list, list_fig_heatmap)

## Subsequent cases / number of imports in each state at every iteration
list_factor_import <- list(fig_005[["factor_import"]], fig_005_wi[["factor_import"]])
## Categories in figure 5
ref_breaks <- c(0, 1, 3, 5, 10, 50)

generate_figure_6(list_factor_import, ref_breaks)

## Figures in supplement
list_out <- list(out_import, out_005, out_09,
                 out_095, out_005_wi, out_095_wi)

supp_fig_param_estimate(list_out, burnin, sample_every)

supp_fig_stratified_state(out_import, dt_cases, burnin, sample_every, max_clust)
supp_fig_stratified_state(out_005, dt_cases, burnin, sample_every, max_clust)
supp_fig_stratified_state(out_095, dt_cases, burnin, sample_every, max_clust)
supp_fig_stratified_state(out_005_wi, dt_cases, burnin, sample_every, max_clust)

list_out <- list(out_import, out_005, out_005_wi)
list_fig <- list(fig_import, fig_005, fig_005_wi, 
                 fig_09, fig_095)


supp_fig_sec_overall(list_out, burnin, sample_every)
supp_fig_sec_map(list_out, dt_cases)
supp_fig_distance_transmission(list_out = list_out, burnin = burnin, 
                               sample_every = sample_every, dt_cases = dt_cases)

supp_fig_eval_cluster(list_fig, titles = NULL)

list_out <- list(out_import, out_005, out_005_wi, 
                 out_09, out_095)
supp_post(list_out, burnin, sample_every)
supp_desc_data(dt_cases = dt_cases)

fig_hist_list <- list(fig_no_gen_sim, fig_005_sim, fig_all_gen_sim, fig_sim)
list_fig_heatmap <- list(fig_no_gen_sim, fig_all_gen_sim)

generate_figure_4_5(fig_hist_list, list_fig_heatmap)

fig_hist_list <- list(fig_no_like, fig_time, fig_time_gen, fig_time_age, fig_time_spa, fig_ref)
list_fig_heatmap <- list(fig_time_gen, fig_time_spa)

generate_figure_4_5(fig_hist_list, list_fig_heatmap)

```

