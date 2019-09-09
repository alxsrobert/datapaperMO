#' Title: Generate the values to describe the inferred clusters and generate the figures
#'
#' @param out: outbreaker_chains data.frame object, output from outbreaker. Contains
#' the inferred transmission chains, parameter estimates, generation numbers and 
#' infection time for every case.
#' 
#' @param dt_cases: Data table. Epi description of the cases. Must contain: the ID, 
#' cluster, State, size_cluster and import status for each case.
#' 
#' @param burnin: Numeric: Length of the burnin period.
#' 
#' @param sample_every: Numeric: Thinning parameter.
#' 
#' @param max_clust: Numeric: Maximum cluster size.
#' 
#' @param thresh_barplot: Vector: Breaks for the barplot
#' 
#' @param diff: Numeric: Length of the heatmap categories.
#' 
#' @return
#' List including \code{med_size_cluster} vector containing the Median values of 
#' cluster size distributions in each run or data included in the histogram ;
#' \code{up_size_cluster} vector containing the 97.5% quantile of cluster size 
#' distributions in each run or data included in the histogram; \code{low_size_cluster}
#' vector containing the 2.5% quantile of cluster size distributions in each run 
#' or data included in the histogram. \code{singletons} vector containing the 
#' 97.5% quantile of cluster size distributions in each run or data included in
#' the histogram. \code{med_imports} vector containing the median number of 
#' imports, and the number of imports correctly inferred in each run or data
#' included in the histogram. \code{up_imports} vector containing the 97.5% 
#' quantile of the number of imports, and the number of imports correctly 
#' inferred in each run or data included in the histogram. \code{low_imports}
#' vector containing the 97.5% quantile of the number of imports, 
#' and the number of imports correctly inferred in each run or data included in the
#' histogram. \code{list_heatmap} vector containing the heatmap data tables for each run.
#' 
#' @export
#'
#' @examples
prepare_for_figures <- function(out, dt_cases, burnin, sample_every, 
                                max_clust, thresh_barplot, diff){
  # Barplot names
  names_barplot <- sapply(1:length(thresh_barplot), function(X){
    if((X + 1) > length(thresh_barplot))
      return(paste0(as.character(thresh_barplot[X]), "+"))
    if(thresh_barplot[X+1] == (thresh_barplot[X] + 1))
      return(as.character(thresh_barplot[X]))
    else 
      return(paste0(as.character(thresh_barplot[X]), "-",
                    as.character(thresh_barplot[X+1] - 1)))
  })
  # Groups_barplot: Correspondence cluser size and group (1 => Group 1; 2 => G2; 3,4 =>G3)
  groups_barplot <- sapply(1:max_clust, function(X){
    return(max(which(thresh_barplot <= X)))
  })
  
  # Breaks for sensitivity and precision vectors
  sensitivity <- seq(diff/2, 1 - diff/2, diff)
  precision <- sensitivity
  
  ## Prepare histogram ##
  # clust_matrix: Ancestor of the tree each case belong to (Column name = infected case, 
  # value = ancestor)
  clust_matrix <- t(apply(out[(burnin/sample_every):dim(out)[1],
                                  grep("alpha", colnames(out))], 1, 
                          function(X){
                            while(any(!is.na(X[X]))){
                              X[!is.na(X[X])] <- X[X[!is.na(X[X])]]
                              
                            }
                            X[!is.na(X)] <- names(X[X[!is.na(X)]])
                            X[is.na(X)] <- names(X[is.na(X)])
                            return(X)
                          }))
  # table_tot: Cluster size distribution
  table_tot <- t(apply(clust_matrix, 1, function(X){
    table_clust <- numeric(max_clust)
    names(table_clust) <- 1:max_clust
    table_clust[names(table(table(X)))] <- table(table(X))
    return(table_clust)
  }))
  # sum_table_tot: Total number of cluster per simulation
  sum_table_tot <- apply(table_tot, 1, sum)
  # Proportion of clusters of each size
  prop_table_tot <- t(apply(table_tot, 1, function(X) return(X)))
  
  
  # dt_singleton: epi description of singletons in reference data
  dt_singletons <- dt_cases[size_cluster == 1,]
  # table_tot: Cluster size distribution considering only singletons
  table_singletons <- t(apply(clust_matrix, 1, function(X){
    table_clust <- numeric(max_clust)
    names(table_clust) <- 1:max_clust
    table_X <- table(X)
    X_singletons <- X[which(is.element(dt_cases$ID,
                                       dt_singletons$ID))]
    table_X_singletons <- table(X)[unique(X_singletons)]
    table_clust[names(table(table_X_singletons))] <- table(table_X_singletons)
    return(table_clust)
  }))
  
  # Group cluster size distributions using groups_barplot
  tot_size_cluster_barplot <- t(rowsum(x = t(prop_table_tot), 
                                       group = groups_barplot))
  colnames(tot_size_cluster_barplot) <- names_barplot
  # Same with only singletons
  size_cluster_singletons <- t(rowsum(x = t(table_singletons), 
                                      group = groups_barplot))
  colnames(tot_size_cluster_barplot) <- names_barplot
  
  # Median and 95% CI for cluster size distribution
  med_prop_singletons <- apply(size_cluster_singletons, 2, median)
  med_size_cluster_barplot <- apply(tot_size_cluster_barplot, 2, median)
  low_size_cluster_barplot <- apply(tot_size_cluster_barplot, 2, function(X) 
    return(quantile(x = X, probs = 0.025)))
  up_size_cluster_barplot <- apply(tot_size_cluster_barplot, 2, function(X) 
    return(quantile(x = X, probs = 0.975)))
  
  ## Imports status ###
  # new: inmput:  ancestor matrix, 
  # Output: number of imports, and number of imports that were imports in the reference
  new <- function(clust_matrix){
    clust_matrix_nb <- gsub(pattern = "alpha_", 
                            replacement = "", 
                            clust_matrix)
    colnames(clust_matrix_nb) <- gsub(pattern = "alpha_", replacement = "", 
                                      colnames(clust_matrix))
    new <- t(apply(clust_matrix_nb, 1, function(X){ 
      # Import at this iteration
      import_iter <- which(X == names(X))
      # Import in data
      index <- length(which(dt_cases[import_iter, import == TRUE]))
      return(c("Imports" = length(import_iter),
               "Import status \ncorrectly inferred" = index
      ))
    }))
    return(new)
  }
  import_infer <- new(clust_matrix = clust_matrix)
  
  # Median and 95% CI for number of imports
  med <- apply(import_infer, 2, median)
  low <- apply(import_infer, 2, function(X) quantile(X, 0.025))
  up <- apply(import_infer, 2, function(X) quantile(X, 0.975))
  
  ## Prepare heatmap ##
  
  clust_matrix_nb <- gsub(pattern = "alpha_", replacement = "", 
                          clust_matrix)
  names(clust_matrix_nb) <- gsub(pattern = "alpha_", replacement = "", 
                                 names(clust_matrix))
  # Which cases are not singletons
  not_singletons <- which(dt_cases$size_cluster>1)
  # Proportion of the reference cluster in the inferred cluster
  prop_clust_in_inferred <- sapply(not_singletons, function(X){
    clust_X <- dt_cases[X, cluster]
    ref_clust <- which(dt_cases$cluster == clust_X)
    clust_matrix_X <- clust_matrix_nb[, ref_clust]
    dim_X <- length(ref_clust)-1
    tab_X <- (apply(clust_matrix_X == clust_matrix_nb[, X], 1, sum)-1)/dim_X
    return(tab_X)})
  
  # Median and 95% CI for sensitivity
  med_clust_in_inferred <- apply(prop_clust_in_inferred, 2, median)
  low_clust_in_inferred <- apply(prop_clust_in_inferred, 2, function(X) 
    quantile(X, 0.025))
  up_clust_in_inferred <- apply(prop_clust_in_inferred, 2, function(X) 
    quantile(X, 0.975))
  
  # Proportion of the inferred cluster in the reference cluster
  prop_inferred_in_clust <- sapply(not_singletons, function(X){
    clust_X <- dt_cases[X, cluster]
    ref_clust <- which(dt_cases$cluster == clust_X)
    clust_matrix_X <- clust_matrix_nb == clust_matrix_nb[, X]
    den <- (apply(clust_matrix_X, 1, sum)-1)
    clust_matrix_X <- clust_matrix_X[,-ref_clust]
    misplaced <- apply(clust_matrix_X, 1, sum) / den
    misplaced[den == 0] <- 0
    return(misplaced)
  })
  # Median and 95% CI for precision
  med_inferred_in_clust <- apply(prop_inferred_in_clust, 2, median)
  up_inferred_in_clust <- apply(prop_inferred_in_clust, 2, 
                                function(X) quantile(X, 0.025))
  low_inferred_in_clust <- apply(prop_inferred_in_clust, 2, 
                                 function(X) quantile(X, 0.975))
  
  dt_prop <- cbind.data.frame(case = not_singletons, 
                              sensitivity = med_clust_in_inferred, 
                              precision = 1 - med_inferred_in_clust)
  dt_heatmap <- cbind(sensitivity = rep(sensitivity, length(precision)),
                      precision = rep(precision, each = length(sensitivity)))
  dt_heatmap <- as.data.table(dt_heatmap)
  dt_heatmap$number <- apply(dt_heatmap, 1, function(X){
    min_sens <- X["sensitivity"]-diff/2
    max_sens <- X["sensitivity"]+diff/2
    min_prec <- X["precision"]-diff/2
    max_prec <- X["precision"]+diff/2
    
    if(max_sens == 1 & max_prec == 1){
      return(length(which(dt_prop$sensitivity >= min_sens & 
                            dt_prop$precision >= min_prec)))
    } else{
      if(max_sens == 1){
        return(length(which(dt_prop$sensitivity >= min_sens & 
                              dt_prop$precision >= min_prec & 
                              dt_prop$precision < max_prec)))
      } 
      if(max_prec == 1){
        return(length(which(dt_prop$sensitivity >= min_sens & 
                              dt_prop$sensitivity < max_sens &
                              dt_prop$precision >= min_prec)))
      } 
    }
    return(length(which(dt_prop$sensitivity >= min_sens & 
                          dt_prop$sensitivity < max_sens &
                          dt_prop$precision >= min_prec &
                          dt_prop$precision < max_prec)))
  })
  dt_heatmap$prop <- dt_heatmap$number/sum(dt_heatmap$number)

  # Proportion subsequent cases / imports per state
  factor_import <- case_state(clust_matrix, dt_cases)
  
  
  return(list(low_size_cluster_barplot = low_size_cluster_barplot,
              med_size_cluster_barplot = med_size_cluster_barplot,
              up_size_cluster_barplot = up_size_cluster_barplot,
              
              med_prop_singletons = med_prop_singletons,
              
              low = low, med = med, up = up,
              
              dt_heatmap = dt_heatmap, factor_import = factor_import,
              
              med_sensitivity =  med_clust_in_inferred,
              low_sensitivity = low_clust_in_inferred,
              up_sensitivity = up_clust_in_inferred,
              
              med_precision = 1 - med_inferred_in_clust,
              low_precision = 1 - low_inferred_in_clust,
              up_precision = 1 - up_inferred_in_clust))
  
}

#' Title: Generate the number of subsequent cases per import in each state
#'
#' @param clust_matrix: Matrix: ancestor of each case
#' @param dt_cases: Reference dataset
#'
#' @return
#' @param factor_import Proportion subsequent cases / imports per state
#'
#' @export
#'
#' @examples
case_state <- function(clust_matrix, dt_cases){
  clust_matrix_nb <- gsub(pattern = "alpha_", 
                          replacement = "", 
                          clust_matrix)
  colnames(clust_matrix_nb) <- gsub(pattern = "alpha_", replacement = "", 
                                    colnames(clust_matrix))
  factor_import <- t(apply(clust_matrix_nb, 1, function(X){
    import_iter <- which(X == names(X))
    factor_import <- table(dt_cases[!import_iter, State])/
      table(dt_cases[import_iter, State])
    factor_import[!is.finite(factor_import)] <- NA
    return(factor_import)
  }))
  return(factor_import)
}


