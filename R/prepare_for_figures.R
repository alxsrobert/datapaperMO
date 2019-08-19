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


prepare_for_figures <- function(out, dt_cases, burnin, sample_every, max_clust,
                                thresh_barplot, diff){
  
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
  
  
  sensitivity <- seq(diff/2, 1 - diff/2, diff)
  precision <- ref
  
  clust_matrix <- t(apply(out[(burnin/sample_every):dim(out)[1],
                                  grep("alpha", colnames(out))], 1, function(X){
                                    while(any(!is.na(X[X]))){
                                      X[!is.na(X[X])] <- X[X[!is.na(X[X])]]
                                      
                                    }
                                    X[!is.na(X)] <- names(X[X[!is.na(X)]])
                                    X[is.na(X)] <- names(X[is.na(X)])
                                    return(X)
                                  }))
  table_tot <- t(apply(clust_matrix, 1, function(X){
    table_clust <- numeric(max_clust)
    names(table_clust) <- 1:max_clust
    table_clust[names(table(table(X)))] <- table(table(X))
    return(table_clust)
  }))
  
  
  ## Prepare histogram ##
  
  dt_singletons <- dt_cases[size_cluster == 1,]
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
  
  sum_table_tot <- apply(table_tot, 1, sum)
  prop_table_tot <- t(apply(table_tot, 1, function(X) return(X)))
  
  tot_size_cluster_barplot <- t(rowsum(x = t(prop_table_tot), 
                                       group = groups_barplot))
  colnames(tot_size_cluster_barplot) <- names_barplot
  
  size_cluster_singletons <- t(rowsum(x = t(table_singletons), 
                                      group = groups_barplot))
  colnames(tot_size_cluster_barplot) <- names_barplot
  
  med_prop_singletons <- apply(size_cluster_singletons, 2, median)
  med_size_cluster_barplot <- apply(tot_size_cluster_barplot, 2, median)
  low_size_cluster_barplot <- apply(tot_size_cluster_barplot, 2, function(X) 
    return(quantile(x = X, probs = 0.025)))
  up_size_cluster_barplot <- apply(tot_size_cluster_barplot, 2, function(X) 
    return(quantile(x = X, probs = 0.975)))
  
  ## Imports status ###
  
  new <- function(clust_matrix){
    clust_matrix_nb <- gsub(pattern = "alpha_", 
                            replacement = "", 
                            clust_matrix)
    colnames(clust_matrix_nb) <- gsub(pattern = "alpha_", replacement = "", 
                                      colnames(clust_matrix))
    new <- t(apply(clust_matrix_nb, 1, function(X){
      import_iter <- which(X == names(X))
      index <- length(which(dt_cases[import_iter, import == TRUE]))
      not_import <- length(which(dt_cases[import_iter, import == FALSE] &
                                   dt_cases[import_iter, import == 2]))
      missed <- length(which(!is.element(which(dt_cases$import == 1),
                                         import_iter)))
      return(c("Imports" = length(import_iter),
               "Import status \ncorrectly inferred" = index
      ))
    }))
    return(new)
  }
  import_infer <- new(clust_matrix = clust_matrix)
  
  med <- apply(import_infer, 2, median)
  low <- apply(import_infer, 2, function(X) quantile(X, 0.025))
  up <- apply(import_infer, 2, function(X) quantile(X, 0.975))
  
  ## Prepare heatmap ##
  
  clust_matrix_nb <- gsub(pattern = "alpha_", replacement = "", 
                          clust_matrix)
  names(clust_matrix_nb) <- gsub(pattern = "alpha_", replacement = "", 
                                 names(clust_matrix))
  
  not_singletons <- which(dt_cases$size_cluster>1)
  prop_clust_in_inferred <- sapply(not_singletons, function(X){
    clust_X <- dt_cases[X, cluster]
    ref_clust <- which(dt_cases$cluster == clust_X)
    clust_matrix_X <- clust_matrix_nb[, ref_clust]
    dim_X <- length(ref_clust)-1
    tab_X <- (apply(clust_matrix_X == clust_matrix_nb[, X], 1, sum)-1)/dim_X
    return(tab_X)})
  
  med_clust_in_inferred <- apply(prop_clust_in_inferred, 2, median)
  low_clust_in_inferred <- apply(prop_clust_in_inferred, 2, function(X) 
    quantile(X, 0.025))
  up_clust_in_inferred <- apply(prop_clust_in_inferred, 2, function(X) 
    quantile(X, 0.975))
  
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
    if(X["sensitivity"] == 1-diff/2){
      return(length(which(dt_prop$sensitivity>=X["sensitivity"]-diff/2 & 
                            dt_prop$sensitivity<=X["sensitivity"]+diff/2 &
                            dt_prop$precision>=X["precision"]-diff/2 & 
                            dt_prop$precision<X["precision"]+diff/2)))    
    }
    if(X["precision"] == 1-diff/2){
      return(length(which(dt_prop$sensitivity>=X["sensitivity"]-diff/2 & 
                            dt_prop$sensitivity<X["sensitivity"]+diff/2 &
                            dt_prop$precision>=X["precision"]-diff/2 & 
                            dt_prop$precision<=X["precision"]+diff/2)))    
    }
    return(length(which(dt_prop$sensitivity>=X["sensitivity"]-diff/2 & 
                          dt_prop$sensitivity<X["sensitivity"]+diff/2 &
                          dt_prop$precision>=X["precision"]-diff/2 &
                          dt_prop$precision<X["precision"]+diff/2)))
  })
  dt_heatmap$prop <- dt_heatmap$number/sum(dt_heatmap$number)
  # y_max <- max(dt_heatmap[prop>0, precision])
  # x_max <- max(dt_heatmap[prop>0, sensitivity])
  # dt_heatmap <- dt_heatmap[precision <= y_max & sensitivity <= x_max,]

  factor_import <- case_state(clust_matrix, dt_cases)
  
  
  clust_transmission <- t(apply(out[(burnin/sample_every):dim(out)[1],
                                    grep("alpha", colnames(out))], 1, function(X){
                                      X[!is.na(X)] <- names(X[X[!is.na(X)]])
                                      X[is.na(X)] <- names(X[is.na(X)])
                                      return(X)
                                    }))
  transmission_distance <- t(apply(clust_transmission, 1, function(X){
    X <- X[X != names(X)]
    X <- gsub(pattern = "alpha_", replacement = "", X)
    names(X) <- gsub(pattern = "alpha_", replacement = "", names(X))
    dist <- dt_distance[paste0(dt_us_cases[as.numeric(X), INCITS], 
                               "_", dt_us_cases[as.numeric(names(X)), INCITS]),dist]
    h <- hist(dist, breaks = c(0, 10, 20, 50, 100), plot = F)
    trans_dist <- h$counts
    names(trans_dist) <- c("0-10", "10-20", "20-50", "50-100")
    return(trans_dist/sum(trans_dist))
  }))
  med_dist <- apply(transmission_distance, 2, median)
  low_dist <- apply(transmission_distance, 2, function(X) 
    return(quantile(x = X, probs = 0.025)))
  up_dist <- apply(transmission_distance, 2, function(X) 
    return(quantile(x = X, probs = 0.975)))
  
  
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
              up_precision = 1 - up_inferred_in_clust, 
              
              med_dist = med_dist, low_dist = low_dist, up_dist = up_dist))
  
}
