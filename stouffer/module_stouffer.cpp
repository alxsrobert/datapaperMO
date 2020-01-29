// [[Rcpp::depends(measlesoutbreaker)]]

#include <Rcpp.h>
#include <Rmath.h>
#include "module_stouffer.h"
#include "measlesoutbreaker.h"

// IMPORTANT: ON INDEXING VECTORS AND ANCESTRIES

// Most of the functions implemented here are susceptible to be called from R
// via Rcpp, and are therefore treated as interfaces. This causes a number of
// headaches when using indices of cases defined in R (1:N) to refer to elements
// in Rcpp / Cpp vectors (0:N-1). By convention, we store all data on the
// original scale (1:N), and modify indices whenever accessing elements of
// vectors. In other words, in an expression like 'alpha[j]', 'j' should always
// be on the internal scale (0:N-1).

// In all these functions, 'SEXP i' is an optional vector of case indices, on
// the 1:N scale.



//  This function compute the spatial log likelihood distribution from parameters a 
//  and b

Rcpp::List cpp_stouffer_matrix(Rcpp::NumericVector population, Rcpp::NumericMatrix distance,
                               Rcpp::NumericMatrix ances, double a, 
                               double gamma, int nb_cases) {
  
  int size_pop = population.size();
  Rcpp::NumericVector population_a(size_pop);
  // log_s_dens with no missing generation
  Rcpp::NumericMatrix probs(nb_cases, nb_cases);
  // log_s_dens with one missing generation
  Rcpp::NumericMatrix probs2(nb_cases, nb_cases);
  Rcpp::NumericVector sum_pop(size_pop);
  Rcpp::NumericMatrix nb_move(size_pop, size_pop);
  double thresh_probs = 0.0;
  if(size_pop <1000) thresh_probs = 0.000001;
  else if(size_pop <3000) thresh_probs = 0.00001;
  
  int j, k, l;
  for(k = 0; k<size_pop; k++){
    population_a[k] = pow(population[k], a);
  }

  for(k = 0; k < size_pop; k++){
    for(j = 0; j < size_pop; j++){
      if(distance(k,j) < gamma){
        nb_move(k,j) = population_a[k]/pow(distance(k,j), a);
        sum_pop[j] += nb_move(k,j);
        if(k != j){
          nb_move(j,k) = nb_move(k,j) * population_a[j]/population_a[k];
          sum_pop[k] += nb_move(j,k);
        }
      }
    }
    for(j = 0; j<size_pop; j++){
      if(distance(k,j) <= gamma){
        if(k < nb_cases && j < nb_cases){
          probs(j, k) = nb_move(j, k)/sum_pop[k];
        }
        else nb_move(j, k) = nb_move(j, k) / sum_pop[k];
      }
    }
  }
  for(k = 0; k<nb_cases; k++){
    for(j = 0; j<nb_cases; j++){
      if(ances(k, j) == 1 && distance(k,j) <= gamma){
        for(l = 0; l < size_pop; l++){
          if(l < nb_cases){
            if((probs(k, l) * probs(l, j)) > thresh_probs &&
               distance(k, l) <= gamma && distance(l, j) <= gamma)
              probs2(k,j) += probs(k, l) * probs(l, j);
          }
          else
            if((nb_move(k, l) * nb_move(l, j)) > thresh_probs &&
               distance(k, l) <= gamma && distance(l, j) <= gamma)
              probs2(k,j)  += nb_move(k, l) * nb_move(l, j);
        }
      }
    }
  }
  for(k = 0; k<nb_cases; k++){
    for(j = 0; j<nb_cases; j++){
      if(ances(k, j) == 1 && distance(k,j) <= gamma){
        probs2(j, k) = log(probs2(j, k));
        probs(j, k) = log(probs(j, k));
      } else if(distance(k,j) > gamma){
        probs2(j, k) = -1000;
        probs(j, k) = -1000;
      }
    }
  }
  
  Rcpp::List new_log_s_dens = Rcpp::List::create(probs, probs2);
  
  return(new_log_s_dens);
}

Rcpp::List cpp_stouffer_move_a(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                               Rcpp::RObject custom_ll, Rcpp::RObject custom_prior) {
  // Import
  Rcpp::List new_param = clone(param);
  double gamma = config["gamma"];
  
  Rcpp::IntegerVector region = data["region"];
  Rcpp::NumericMatrix distance = data["distance"];
  Rcpp::NumericVector population = data["population"];
  Rcpp::NumericVector limits = config["prior_a"];
  Rcpp::NumericMatrix can_be_ances_reg = data["can_be_ances_reg"];
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
  // printf("YO1\n");
  new_param["log_s_dens"] = cpp_stouffer_matrix(population, distance, can_be_ances_reg,
                                      new_a[0], gamma, nb_cases);
  // printf("YO2\n");
  // compute likelihoods
  old_logpost = measlesoutbreaker::cpp_ll_space(data, config, param, R_NilValue, custom_ll);
  new_logpost = measlesoutbreaker::cpp_ll_space(data, config, new_param, R_NilValue, custom_ll);
  
  // compute priors
  
  old_logpost += measlesoutbreaker::cpp_prior_a(param, config, custom_prior);
  new_logpost += measlesoutbreaker::cpp_prior_a(new_param, config, custom_prior);
  // acceptance term
  
  p_accept = exp(new_logpost - old_logpost);
  
  
  // acceptance: the new value is already in a, so we only act if the move is
  // rejected, in which case we restore the previous ('old') value
  
  if (p_accept < unif_rand()) { // reject new values
    return param;
  }
  
  return new_param;
}

