prepare_for_figures <- function(out, dt_cases){
  
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
  
  tot_size_cluster_barplot <- t(rowsum(x = t(prop_table_tot), group = groups_barplot))
  colnames(tot_size_cluster_barplot) <- names_barplot
  
  size_cluster_singletons <- t(rowsum(x = t(table_singletons), group = groups_barplot))
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
      # import <- length(which(dt_cases[import_iter, import == FALSE] &
      #                          dt_cases[import_iter, import == 1]))
      # index <- index + import
      not_import <- length(which(dt_cases[import_iter, import == FALSE] &
                                   dt_cases[import_iter, import == 2]))
      missed <- length(which(!is.element(which(dt_cases$import == 1),
                                         import_iter)))# / length(which(dt_cases$import == 1))
      return(c("Imports" = length(import_iter),
               "Import status \ncorrectly inferred" = index#,#/length(import_iter),
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
  up_inferred_in_clust <- apply(prop_inferred_in_clust, 2, function(X) quantile(X, 0.025))
  low_inferred_in_clust <- apply(prop_inferred_in_clust, 2, function(X) quantile(X, 0.975))
  
  dt_prop <- cbind.data.frame(case = not_singletons, 
                              ref = med_clust_in_inferred, 
                              missed = med_inferred_in_clust)
  dt_heatmap <- cbind(ref = rep(ref, length(missed)),
                      missed = rep(missed, each = length(ref)))
  dt_heatmap <- as.data.table(dt_heatmap)
  dt_heatmap$number <- apply(dt_heatmap, 1, function(X){
    if(X["ref"] == 1-diff/2){
      return(length(which(dt_prop$ref>=X["ref"]-diff/2 & 
                            dt_prop$ref<=X["ref"]+diff/2 &
                            dt_prop$missed>=X["missed"]-diff/2 & 
                            dt_prop$missed<X["missed"]+diff/2)))    
    }
    if(X["missed"] == 1-diff/2){
      return(length(which(dt_prop$ref>=X["ref"]-diff/2 & 
                            dt_prop$ref<X["ref"]+diff/2 &
                            dt_prop$missed>=X["missed"]-diff/2 & 
                            dt_prop$missed<=X["missed"]+diff/2)))    
    }
    return(length(which(dt_prop$ref>=X["ref"]-diff/2 & 
                          dt_prop$ref<X["ref"]+diff/2 &
                          dt_prop$missed>=X["missed"]-diff/2 &
                          dt_prop$missed<X["missed"]+diff/2)))
  })
  dt_heatmap$prop <- dt_heatmap$number/sum(dt_heatmap$number)
  y_max <- max(dt_heatmap[prop>0, missed])
  x_max <- max(dt_heatmap[prop>0, ref])
  dt_heatmap <- dt_heatmap[missed <= y_max & ref <= x_max,]
  
  return(list(low_size_cluster_barplot = low_size_cluster_barplot,
              med_size_cluster_barplot = med_size_cluster_barplot,
              up_size_cluster_barplot = up_size_cluster_barplot,
              med_prop_singletons = med_prop_singletons,
              low = low,
              med = med,
              up = up,
              dt_heatmap = dt_heatmap))
  
}