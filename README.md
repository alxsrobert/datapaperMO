# datapaperMO
analysis of toy_outbreak dataset from the measlesoutbreaker package

**install package measlesoutbreaker from "alxsrobert/measlesoutbreaker". To generate all figures from main and supplement, run the script generate_figure.R.**

```R
install.packages("devtools")
library(devtools)
install_github("alxsrobert/measlesoutbreaker")
library(measlesoutbreaker)
source("R/generate_figures.R")
```

**To generate a new simulated dataset with different input parameters (number of cases, R0, spatial parameters...), edit the lines 10-37 in analysis_generated_data.R. Run the rest of the script to generate the measlesoutbreaker runs on the dataset.**

**List of files in the repository:**
* In data: 
  * State_initials.csv	: Data table with state names and acronyms.
  * pop_center.csv: Population centroid and population of every county (Source: https://www.census.gov/geographies/reference-files/2010/geo/2010-centers-population.html (per state : UNITED STATES)).

* In R:
  * library_importation.R: import libraries for analysis.
  * function_generate_dataset.R: Functions to generate toy_outbreak.
  * function_prepare_for_figures.R: Function to generate summary statistics on run.
  * function_generate_figures_main.R: function to generate figures similar to figure 3, 4, 5 and 6.
  * function_supplement_figures.R: function to generate figures similar to figures in the supplement.
  * load_analysis_threshold.R: Load measlesoutbreaker runs generated with different import thresholds.
  * load_analysis_likelihoods.R: Load measlesoutbreaker runs generated with different likelihoods.
  * load_analysis_genotype.R: Load measlesoutbreaker runs generated with different proportion of genotyped cases.
  * load_data_distributions.R: Load prior distributions.
  * analysis_generated_data.R: Script to run measlesoutbreaker on toy_outbreak.
  * generate_figure.R: Generate all figures (main + supplement).
