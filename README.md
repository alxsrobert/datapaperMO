# datapaperMO
analysis of fake_outbreak dataset from the measlesoutbreaker package

In data: 
- State_initials.csv	: Data table with state names and acronyms.
- pop_center.csv: Population centroid and population of every county (Source: https://www.census.gov/geographies/reference-files/2010/geo/2010-centers-population.html (per state : UNITED STATES)).
- social_contact_uk.csv : Contact Matrix of contacts between age groups (Source: Social Contacts and Mixing Patterns Relevant to the Spread of Infectious Diseases)

In src:
- library_importation.R: import libraries for analysis
- load_data_distributions.R: Get prior distributions
- generate_dataset.R: Functions to generate toy_outbreak
- prepare_for_figures.R: Function to generate summary statistics on run.
- function_generate_figures_main.R: function to generate figures similar to figure 3, 4 and 5
- analysis_generated_data.R: Script to run measlesoutbreaker on toy_outbreak  
- load_analysis_for_figure.R: Script to analyse the output of measlesoutbreaker runs
