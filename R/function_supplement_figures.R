supp_fig_eval_cluster <- function(list_fig){
  par(mfrow = c(length(list_fig),2))
  lapply(list_fig, function(X){
    low_sens <- X$low_sensitivity[order(X$med_sensitivity)]
    up_sens <- X$up_sensitivity[order(X$med_sensitivity)]
    med_sens <- X$med_sensitivity[order(X$med_sensitivity)]
    
    low_prec <- X$low_precision[order(X$med_precision)]
    up_prec <- X$up_precision[order(X$med_precision)]
    med_prec <- X$med_precision[order(X$med_precision)]
    
    
    plot(med_sens, type = "l", ylab = "", xlab = "")
    xx=c(1:length(med_sens), 
         rev(1:length(med_sens)))
    yy=c(low_sens, 
         rev(up_sens))
    polygon(xx, yy, col = transp("black", 0.3), border = NA)
    
    plot(med_prec, type = "l", ylab = "", xlab = "")
    xx=c(1:length(med_prec), 
         rev(1:length(med_prec)))
    yy=c(low_prec, 
         rev(up_prec))
    polygon(xx, yy, col = transp("black", 0.3), border = NA)
  })
}

supp_estim_params <- function(out){
  pi = out$pi
  med_pi <- median(pi)
  up_pi <- quantile(pi, 0.975) %>% as.numeric
  low_pi <- quantile(pi, 0.025) %>% as.numeric
  
  a = out$a
  med_a <- median(a)
  up_a <- quantile(a, 0.975) %>% as.numeric
  low_a <- quantile(a, 0.025) %>% as.numeric
  
  b = out$b
  med_b <- median(b)
  up_b <- quantile(b, 0.975) %>% as.numeric
  low_b <- quantile(b, 0.025) %>% as.numeric
  
  return(c(pi = med_pi, pi_low = low_pi, pi_up = up_pi,
           a = med_a, a_low = low_a, a_up = up_a,
           b = med_b, b_low = low_b, b_up = up_b))
}

supp_fig_param_estimate <- function(list_out){
  output <- do.call(rbind,lapply(list_out, estim_params))
  
  par(fig = c(0, 1, 0.68, 1), family = "sans", omi = c(0,0,0,0), cex.lab = 1.4,
      cex.axis = 1.2, mar = c(0,4,-0.5,1)+ 1.1, las=1, bty="l", tcl=-1, lwd=2,
      mgp = c(3.5, 1, 0))
  
  plot(output[,"pi"], pch = 17, ylim = c(min(output[,"pi_low"]), 
                                         max(output[, "pi_up"])), 
       ylab  = "Report ratio", xaxt = "n", xlab = "")
  arrows(x0 = 1:length(list_out), y0 = output[,"pi_low"], 
         x1 = , y1 = output[,"pi_up"], angle = 90, 
         code = 3,length = 0.1, lwd =  0.5)
  axis(side = 1, at = 1:dim(output)[1], labels = rep("", dim(output)[1]))
  
  par(fig = c(0, 1, 0.36, 0.68), family = "sans", omi = c(0,0,0,0), cex.lab = 1.4,
      cex.axis = 1.2, mar = c(0,4,-0.5,1)+ 1.1, las=1, bty="l", tcl=-1, lwd=2,
      mgp = c(3.5, 1, 0), new = T)
  
  plot(output[,"a"], pch = 17, ylim = c(0.95*min(output[,"a_low"]), 
                                        1.05*max(output[, "a_up"])),
       ylab = "Spatial parameter (a)", xaxt = "n", xlab = "")
  arrows(x0 = 1:length(list_out), y0 = output[,"a_low"], 
         x1 = , y1 = output[,"a_up"], angle = 90, 
         code = 3,length = 0.1, lwd =  0.5)
  axis(side = 1, at = 1:dim(output)[1], labels = rep("", dim(output)[1]))
  
  par(fig = c(0,1, 0, 0.36), family = "sans", omi = c(0,0,0,0), cex.lab = 1.4,
      cex.axis = 1.2, mar = c(3,4,-1.1,1)+ 1.1, las=1, bty="l", tcl=-1, lwd=2,
      mgp = c(3.5, 1, 0), new = T)
  
  plot(output[,"b"], pch = 17, ylim = c(0.95*min(output[,"b_low"]), 
                                        1.05*max(output[, "b_up"])),
       ylab = "Spatial parameter (b)", xaxt = "n", xlab = "")
  arrows(x0 = 1:length(list_out), y0 = output[,"b_low"], 
         x1 = , y1 = output[,"b_up"], angle = 90, 
         code = 3,length = 0.1, lwd =  0.5)
  
  axis(side = 1, at = 1:dim(output)[1], labels = rep("", dim(output)[1]))
  mgp.axis(side = 1, at = 1:dim(output)[1], labels = rownames(output), mgp = c(3.5,2.5,0))
  
}

supp_fig_stratified_state <- function(out, dt_cases, burnin, sample_every, max_clust){
  states <- names(sort(table(dt_cases$State),decreasing = T))[1:9]

  clust_matrix <- t(apply(out[(burnin/sample_every):dim(out)[1],
                              grep("alpha", colnames(out))], 1, function(X){
                                while(any(!is.na(X[X]))){
                                  X[!is.na(X[X])] <- X[X[!is.na(X[X])]]
                                  
                                }
                                X[!is.na(X)] <- names(X[X[!is.na(X)]])
                                X[is.na(X)] <- names(X[is.na(X)])
                                return(X)
                              }))
  
  par(mfrow = c(3, 3), omi = c(0.2,0.25,0,0), cex.lab = 1.4, cex.axis = 1.2, 
      mar = c(2,1,0,-1)+ 1.1, las=1, bty="l", tcl=-1, lwd=2, mgp = c(3, 1, 0))
  
  for(X in states){
    cases_clust <- which(dt_cases$State == X)
    clust_matrix_clust <- clust_matrix[, cases_clust]

    table_clust <- t(apply(clust_matrix_clust, 1, function(X){
      table_clust <- numeric(max_clust)
      names(table_clust) <- 1:max_clust
      table_clust[names(table(table(X)))] <- table(table(X))
      return(table_clust)
    }))
    size_cluster_inferred_barplot <- 
      cbind("1" = table_clust[,1], "2" = table_clust[,2],
            "3-4" = apply(table_clust[,3:4], 1, sum),
            "5-9" = apply(table_clust[,5:9], 1, sum),
            "10-19" = apply(table_clust[,10:19], 1, sum),
            "20-29" = apply(table_clust[,20:29], 1, sum),
            "30-49" = apply(table_clust[,30:49], 1, sum),
            "50-99" = apply(table_clust[,50:99], 1, sum),
            "100-199" = apply(table_clust[,100:199], 1, sum),
            "200+" = apply(table_clust[,200:dim(table_clust)[2]], 1, sum)
      )
    med_size_cluster_inferred_barplot <- apply(size_cluster_inferred_barplot, 2, median)
    low_size_cluster_inferred_barplot <- apply(size_cluster_inferred_barplot, 2, function(X) 
      return(quantile(x = X, probs = 0.025)))
    up_size_cluster_inferred_barplot <- apply(size_cluster_inferred_barplot, 2, function(X) 
      return(quantile(x = X, probs = 0.975)))
    
    ref_size_clusters <- dt_cases[!duplicated(cluster) & cluster != "." & 
                                    State == X, size_cluster]
    
    ref_size_clusters <- c(rep(1, dim(dt_cases[State == X & cluster == ".",])[1]),
                           ref_size_clusters)
    ref_size_clusters <- factor(ref_size_clusters, levels = 1:max_clust)
    ref_size_clusters <- table(ref_size_clusters)
    
    size_cluster_ref_barplot <- 
      c("1" = as.numeric(ref_size_clusters[1]),
        "2" = as.numeric(ref_size_clusters[2]),
        "3-4" = sum(ref_size_clusters[3:4]),
        "5-9" = sum(ref_size_clusters[5:9]),
        "10-19" = sum(ref_size_clusters[10:19]),
        "20-29" = sum(ref_size_clusters[20:29]),
        "30-49" = sum(ref_size_clusters[30:49]),
        "50-99" = sum(ref_size_clusters[50:99]),
        "100-199" = sum(ref_size_clusters[100:199]),
        "200+" = sum(ref_size_clusters[200:length(ref_size_clusters)]))
    
    
    b <- barplot(rbind(med_size_cluster_inferred_barplot, 
                       size_cluster_ref_barplot),
                 col = grey.colors(5)[c(2, 5)],
                 beside = T, ylim = c(0, max(c(up_size_cluster_inferred_barplot,
                                               size_cluster_ref_barplot))),
                 main = X, border = NA, cex.main = 2)
    arrows(x0 = b[1, ], y0 = low_size_cluster_inferred_barplot, 
           x1 = , y1 = up_size_cluster_inferred_barplot, angle = 90, 
           code = 3,length = 0.1)
  }
  
  
}

supp_fig_sec_overall <- function(list_out, burnin, sample_every){
  breaks <- c(0, 1, 3, 5, 100)
  names_breaks <- c("0", "1-2", "3-4", "5+")
  
  table_secondary_list <- lapply(list_out, function(out_X){
    transmission_matrix <- t(apply(out_X[(burnin/sample_every):dim(out_X)[1],
                                      grep("alpha", colnames(out_X))], 1, 
                                  function(X){
                                    X[!is.na(X)] <- names(X[X[!is.na(X)]])
                                    return(X)
                                  }))
    max_sec <- max(apply(transmission_matrix, 1, function(X) return(max(table(X)))))
    table_secondary <- t(apply(transmission_matrix, 1, function(X){
      table_sec <- numeric(max_sec+1)
      names(table_sec) <- 0:max_sec
      table_sec[names(table(table(X)))] <- table(table(X))
      table_sec[1] <- length(X)-sum(table_sec)
      table_sec <- table_sec/sum(table_sec)

      table_aggreg <- table(cut(x = table(X), 
                                breaks=breaks,
                                include.lowest = T,right = F,
                                labels=names_breaks))
      table_aggreg[1] <- length(X)-sum(table_aggreg)
      table_aggreg <- table_aggreg/sum(table_aggreg)
      return(table_aggreg)
      
    }))
    return(table_secondary)
    
  })
  
  med_sec_list <- lapply(table_secondary_list, function(X) 
    return(as.data.frame(t(apply(X, 2, median)))))
  med_sec <- bind_rows(med_sec_list) %>% as.matrix
  
  low_sec_list <- lapply(table_secondary_list, function(X) 
    return(as.data.frame(t(apply(X, 2, function(Y) 
      quantile(x = Y, probs = 0.025))))))
  low_sec <- bind_rows(low_sec_list) %>% as.matrix
  
  up_sec_list <- lapply(table_secondary_list, function(X) 
    return(as.data.frame(t(apply(X, 2, function(Y) 
      quantile(x = Y, probs = 0.975))))))
  up_sec <- bind_rows(up_sec_list) %>% as.matrix
  
  
  par(mfrow = c(1, 1), cex.lab = 1.5,
      oma = c(0,1,0,0), cex.axis=1.2, cex.main=1.5, #b l t r
      mar = c(3,3,0,1) + 1.3, las=1, bty="l", tcl=-1, lwd=2,
      mgp = c(3, 1.2, 0))
  b <- barplot(med_sec, 
               col = grey.colors(nrow(med_sec)),
               beside = T, 
               ylim = c(0, max(up_sec)),
               xlab = "Number of secondary cases", ylab = "Proportion", 
               border = NA)
  arrows(x0 = b, y0 = low_sec, x1 = , y1 = up_sec, angle = 90, code = 3,
         length = 0.1)

}

supp_fig_sec_map <- function(list_out, dt_cases){
  breaks_map <- c(0, 0.2, 0.4, 0.6, 0.8, 1)
  categ_map <- c("0-0.2",
                 "0.2-0.4",
                 "0.4-0.6",
                 "0.6-0.8",
                 "0.8-1")
  data(state, package = "datasets")
  dt_map_r0 <- data.table(abb = c(state.abb, "COL"))
  dt_map_r0[, region := tolower(c(state.name, "district of columbia"))]
  setkey(dt_map_r0, abb)
  dt_map_r0[, r0 := 0]
  dt_map_r0 <- dt_map_r0[region != "alaska",]
  dt_map_r0 <- dt_map_r0[region != "hawaii",]
  
  dt_map_r0_list <- lapply(list_out, function(out_X){
    transmission_matrix <- t(apply(out_X[(burnin/sample_every):dim(out_X)[1],
                                         grep("alpha", colnames(out_X))], 1, 
                                   function(X){
                                     X[!is.na(X)] <- names(X[X[!is.na(X)]])
                                     return(X)
                                   }))
    sec_per_state <- t(apply(transmission_matrix, 1, function(X){
      X <- factor(X, levels = names(X))
      tab_X <- table(X)
      dt_test <- dt_cases
      dt_test[, sec := tab_X]
      dt_test <- dt_test[, .(State, sec)]
      dt_test <- dt_test[, lapply(.SD, sum), by = State]
      nb_sec <- dt_test$sec
      names(nb_sec) <- dt_test$State
      nb_sec <- nb_sec[names(table(dt_us_cases$State))]/table(dt_us_cases$State)
      
      return(nb_sec)
    }))
    med_r0_X <- apply(sec_per_state, 2, median)
    
    dt_map_r0_X = as.data.table(dt_map_r0)
    dt_map_r0_X[names(med_r0_X), r0 := med_r0_X]
    dt_map_r0_X <- dt_map_r0_X %>% 
      mutate(category=cut(r0, breaks=breaks_map,include.lowest = T,
                          labels=categ_map))
    
    return(dt_map_r0_X)
  })
  
  dt_map <- bind_rows(dt_map_r0_list, .id = "id")
  all_states <- map_data("state")
  
  Total <- merge(all_states, dt_map, by="region")
  
  p <- ggplot(Total, aes(x=Total$long, y=Total$lat)) +   
    facet_grid(id~.)
  p <- p + geom_polygon(data=Total, aes(group = group
                                        # , fill=as.numeric(Total$fact)
                                        , fill = Total$category),
                        color = "black") +
    scale_fill_manual(values = brewer.pal(n = 5, name = 'Purples'),
                      na.value = "grey"
                      , breaks = categ_map
                      # , guide="colorbar"
    )
  p <- p + theme_bw()  + 
    labs(#fill = "Number of measles r0\nreported per state",
      fill = "",
      x="", y="")
  p <- p + scale_y_continuous(breaks=c()) + scale_x_continuous(breaks=c()) + 
    theme(panel.border =  element_blank(), legend.text=element_text(size=25), 
          strip.text.y = element_blank(),
          legend.title = element_text(size = 18), title = element_text(size = 20))
  p
}

supp_fig_distance_transmission <- function(list_out, burnin, sample_every, dt_distance){
  dt_distance[, id_dist := paste0(county1, "_", county2)]
  setkey(dt_distance, id_dist)
  
  list_transmission_dist <- lapply(list_out, function(out_X){
    clust_transmission <- t(apply(out_X[(burnin/sample_every):dim(out_X)[1],
                                      grep("alpha", colnames(out_X))], 1, function(X){
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
    return(transmission_distance)
  })
  med_dist_list <- lapply(list_transmission_dist, function(X) 
    return(as.data.frame(t(apply(X, 2, median)))))
  med_dist <- bind_rows(med_dist_list) %>% as.matrix
  
  low_dist_list <- lapply(list_transmission_dist, function(X) 
    return(as.data.frame(t(apply(X, 2, function(Y) 
      quantile(x = Y, probs = 0.025))))))
  low_dist <- bind_rows(low_dist_list) %>% as.matrix
  
  up_dist_list <- lapply(list_transmission_dist, function(X) 
    return(as.data.frame(t(apply(X, 2, function(Y) 
      quantile(x = Y, probs = 0.975))))))
  up_dist <- bind_rows(up_dist_list) %>% as.matrix
  
  par(mfrow = c(1,1), family = "sans", omi = c(0,0,0,0), cex.lab = 1.4,
      cex.axis = 1.2,  mar = c(3,3.5,0,1)+ 1.1, las=1, bty="l", tcl=-1, lwd=2,
      mgp = c(3, 1, 0))
  
  b <- barplot(rbind(med_dist), 
               col =grey.colors(nrow(med_dist)),
               beside = T, 
               ylim = c(0, 1),
               xlab = "Distance (km)", ylab = "Proportion", border = NA)
  arrows(x0 = b, y0 = low_dist, x1 = , y1 = up_dist, angle = 90, code = 3,
         length = 0.1)

}