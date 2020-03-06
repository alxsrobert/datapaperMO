// [[Rcpp::depends(o2geosocial)]]

#include <Rcpp.h>
#include <Rmath.h>
#include "o2geosocial.h"

// [[Rcpp::export()]]
Rcpp::List cpp_stouffer_move_a(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                               Rcpp::RObject custom_ll, Rcpp::RObject custom_prior) {
  // Import
  Rcpp::List new_param = clone(param);
  double gamma = config["gamma"];
  int max_kappa = config["max_kappa"];
  
  Rcpp::String spatial = config["spatial_method"];
  Rcpp::IntegerVector region = data["region"];
  Rcpp::NumericMatrix distance = data["distance"];
  Rcpp::NumericMatrix can_be_ances_reg = data["can_be_ances_reg"];
  Rcpp::NumericVector population = data["population"];
  Rcpp::NumericVector limits = config["prior_a"];
  
  Rcpp::List new_log_s_dens = new_param["log_s_dens"];
  
  Rcpp::NumericVector new_a = new_param["a"]; // these are just pointers
  Rcpp::NumericMatrix probs = new_log_s_dens[0];
  
  int nb_cases = pow(probs.size(), 0.5);
  
  double sd_a = static_cast<double>(config["sd_a"]);
  
  double old_logpost = 0.0, new_logpost = 0.0, p_accept = 0.0;
  
  // Move new_a
  // proposal (normal distribution with SD: config$sd_a)
  
  new_a[0] += R::rnorm(0.0, sd_a); // new proposed value
  
  if (new_a[0] < limits[0] || new_a[0] > limits[1]) {
    return param;
  }
  new_param["log_s_dens"] = o2geosocial::cpp_log_like(population, distance, can_be_ances_reg,
                                            new_a[0], new_a[0], max_kappa, gamma, spatial, nb_cases);
  
  // compute likelihoods
  old_logpost = o2geosocial::cpp_ll_space(data, config, param, R_NilValue, custom_ll);
  new_logpost = o2geosocial::cpp_ll_space(data, config, new_param, R_NilValue, custom_ll);
  
  
  // compute priors
  
  old_logpost += o2geosocial::cpp_prior_a(param, config, custom_prior);
  new_logpost += o2geosocial::cpp_prior_a(new_param, config, custom_prior);
  // acceptance term
  
  p_accept = exp(new_logpost - old_logpost);
  
  
  // acceptance: the new value is already in a, so we only act if the move is
  // rejected, in which case we restore the previous ('old') value
  
  if (p_accept < unif_rand()) { // reject new values
    return param;
  }
  
  return new_param;
}

