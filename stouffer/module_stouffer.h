#ifndef STOUFFER_MODULE_STOUFFER_H
#define STOUFFER_MODULE_STOUFFER_H

#include <Rcpp.h>
#include <Rmath.h>

// [[Rcpp::export()]]
Rcpp::List cpp_stouffer_matrix(Rcpp::NumericVector population, Rcpp::NumericMatrix distance,
                               Rcpp::NumericMatrix ances, double a, 
                               double gamma, int nb_cases);
// [[Rcpp::export()]]
Rcpp::List cpp_stouffer_move_a(Rcpp::List param, Rcpp::List data, Rcpp::List config,
                               Rcpp::RObject custom_ll = R_NilValue,
                               Rcpp::RObject custom_prior = R_NilValue);

#endif