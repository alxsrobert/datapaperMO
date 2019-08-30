transp <- function(col, alpha=.5){
  res <- apply(col2rgb(col),2, function(c) 
    rgb(c[1]/255, c[2]/255, c[3]/255, alpha))
  return(res)
}

supp_fig_eval_cluster <- function(list_fig, titles){
  par(mfrow = c(length(list_fig),2), oma = c(2,1,0,1), cex.lab = 2, 
      cex.axis = 1.5, mar = c(1,3.5,1,-1)+ 1.1, las=1, bty="l", tcl=-1, lwd=2,
      mgp = c(3, 1, 0))
  for (i in 1:length(list_fig)){
    X <- list_fig[[i]]
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
    mtext(text = titles[i], adj = 0, cex = 1.5, line = 0.5)
    
    plot(med_prec, type = "l", ylab = "", xlab = "")
    xx=c(1:length(med_prec), 
         rev(1:length(med_prec)))
    yy=c(low_prec, 
         rev(up_prec))
    polygon(xx, yy, col = transp("black", 0.3), border = NA)
  }
  title(outer = T, xlab = "Cases", line = 0.5, cex.lab = 2)
  title(outer = T, ylab = "Sensitivity", line = -1, cex.lab = 2)
  par(fig = c(0.5, 1, 0, 1), new = T)
  plot.new()
  title(ylab = "Precision", line = 3, cex.lab = 2)
  
}

estim_params <- function(out, burnin, sample_every){
  pi = out[(burnin/sample_every):dim(out)[1], "pi"]
  med_pi <- median(pi)
  up_pi <- quantile(pi, 0.975) %>% as.numeric
  low_pi <- quantile(pi, 0.025) %>% as.numeric
  
  a = out[(burnin/sample_every):dim(out)[1], "a"]
  med_a <- median(a)
  up_a <- quantile(a, 0.975) %>% as.numeric
  low_a <- quantile(a, 0.025) %>% as.numeric
  
  b = out[(burnin/sample_every):dim(out)[1], "b"]
  med_b <- median(b)
  up_b <- quantile(b, 0.975) %>% as.numeric
  low_b <- quantile(b, 0.025) %>% as.numeric
  
  return(c(pi = med_pi, pi_low = low_pi, pi_up = up_pi,
           a = med_a, a_low = low_a, a_up = up_a,
           b = med_b, b_low = low_b, b_up = up_b))
}

supp_fig_param_estimate <- function(list_out, burnin, sample_every){
  output <- do.call(rbind,lapply(list_out, function(X)
    return(estim_params(X, burnin, sample_every))))

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
  title(xlab = "Cluster size", outer = T, line = .5, cex.lab = 2)
  title(ylab = "Number", outer = T, line = 0.5, cex.lab = 2)
  par(mfrow = c(1, 1),
      omi = c(0.2,0.25,0,0),
      cex.lab = 1.4,
      cex.axis = 1.2, 
      mar = c(2,1,0,-1)+ 1.1,
      las=1, bty="l", tcl=-1, lwd=2,
      mgp = c(3, 1, 0), new = T)
  legend("center",
         fill = c(grey.colors(5)[c(2, 5)], NA), cex = 1.2,
         border = NA, 
         legend = c("Inferred distribution", 
                    "Epidemiological cluster"), bty = "n")
  
  
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
      nb_sec <- nb_sec[names(table(dt_cases$State))]/table(dt_cases$State)
      
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
  
  p <- ggplot(Total, aes(x=Total$long, y=Total$lat)) + facet_grid(id~.)
  p <- p + geom_polygon(data=Total, aes(group = group, fill = Total$category),
                        color = "black") +
    scale_fill_manual(values = brewer.pal(n = 5, name = 'Purples'),
                      na.value = "grey", breaks = categ_map)
  p <- p + theme_bw()  + labs(fill = "", x="", y="")
  p <- p + scale_y_continuous(breaks=c()) + scale_x_continuous(breaks=c()) + 
    theme(panel.border =  element_blank(), legend.text=element_text(size=25), 
          strip.text.y = element_blank(),
          legend.title = element_text(size = 18), title = element_text(size = 20))
  p
}

supp_fig_distance_transmission <- function(list_out, burnin, sample_every, 
                                           dt_cases){
  # Load excel file with population centroid in every US county
  dt_loca_pop_center <- as.data.table(read.csv(file = "data/pop_center.csv",
                                               stringsAsFactors = FALSE, 
                                               colClasses = c("character", "character",
                                                              "character","character",
                                                              "numeric", "numeric",
                                                              "numeric", "character")))
  # dt_state_county includes county ID and the state it belongs to
  dt_state_county <- dt_loca_pop_center[,.(ID_COUNTY, STNAME)]
  setkey(dt_state_county, ID_COUNTY)
  # Create empty data table to story distance between every UK county
  dt_distance <- as.data.table(matrix(0, nrow = nrow(dt_loca_pop_center)**2,
                                      ncol = 9))
  colnames(dt_distance) <- c("county1", "county2", "distance_km",
                             "pop_county1", "long1", "lat1",
                             "pop_county2", "long2", "lat2")
  # 1st column: county1
  dt_distance[, county1 := as.character(county1)]
  dt_distance[, county1 := rep(dt_loca_pop_center$ID_COUNTY, nrow(dt_loca_pop_center))]
  setkey(dt_loca_pop_center, ID_COUNTY)
  # 5th column: long value of county1's population centroid
  dt_distance[, long1 := dt_loca_pop_center[dt_distance$county1, LONGITUDE]]
  # 6th column: lat value of county1's population centroid
  dt_distance[, lat1 := dt_loca_pop_center[dt_distance$county1,LATITUDE]]
  # 2nd column: county2
  dt_distance[, county2 := rep(dt_loca_pop_center$ID_COUNTY,each = nrow(dt_loca_pop_center))]
  # 8th column: long value of county2's population centroid
  dt_distance[, long2 := dt_loca_pop_center[dt_distance$county2,LONGITUDE]]
  # 9th column: lat value of county2's population centroid
  dt_distance[, lat2 := dt_loca_pop_center[dt_distance$county2, LATITUDE]]
  long1 <- dt_distance$long1
  lat1 <- dt_distance$lat1
  long2 <- dt_distance$long2
  lat2 <- dt_distance$lat2
  mat1 <- matrix(c(long1,lat1), ncol = 2)
  mat2 <- matrix(c(long2, lat2), ncol = 2)
  dist <- numeric(dim(dt_distance)[1])
  # Compute distance between every centroid
  dist <- distGeo(mat1, mat2)/1000
  # 3rd column: distance between counties
  dt_distance[, dist := dist]
  # 4th column: population in county1
  dt_distance[, pop_county1 := dt_loca_pop_center[dt_distance$county1,POPULATION]]
  # 7th column: population in county2
  dt_distance[, pop_county2 := dt_loca_pop_center[dt_distance$county2, POPULATION]]
  dt_distance[, id_dist := paste0(county1, "_", county2)]
  
  setkey(dt_distance, id_dist)
  
  list_transmission_dist <- lapply(list_out, function(out_X){
    clust_transmission <- t(apply(out_X[(burnin/sample_every):dim(out_X)[1],
                                        grep("alpha", colnames(out_X))], 1, 
                                  function(X){
                                    X[!is.na(X)] <- names(X[X[!is.na(X)]])
                                    X[is.na(X)] <- names(X[is.na(X)])
                                    return(X)}))
    
    transmission_distance <- t(apply(clust_transmission, 1, function(X){
      X <- X[X != names(X)]
      X <- gsub(pattern = "alpha_", replacement = "", X)
      names(X) <- gsub(pattern = "alpha_", replacement = "", names(X))
      dist <- dt_distance[paste0(dt_cases[as.numeric(X), county], 
                                 "_", dt_cases[as.numeric(names(X)), county]),dist]
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

supp_desc_data <- function(dt_cases){
  data(state, package = "datasets")
  dt_map_clusters <- data.table(abb = c(state.abb, "COL"))
  dt_map_clusters[, region := tolower(c(state.name, "district of columbia"))]
  setkey(dt_map_clusters, abb)
  dt_map_clusters[, clusters := 0]
  tab_clust_state <- table(dt_cases$State, dt_cases$cluster)
  nb_clust_state <- apply(tab_clust_state, 1, function(X){
    if(any(names(X) == ".")){
      if(X["."]>0)
        return(length(which(X>0))+X["."]-1)
      return(length(which(X>0))+X["."])
    }
    else
      return(length(which(X>0)))
  })
  
  dt_map_clusters[names(nb_clust_state),
                  clusters := as.numeric(nb_clust_state)]
  dt_map_clusters <- dt_map_clusters[region != "alaska",]
  dt_map_clusters <- dt_map_clusters[region != "hawaii",]
  all_states <- map_data("state")
  Total <- merge(all_states, dt_map_clusters, by="region")
  
  dt_cases[, age_group_char := paste0((age_group-1)*5, "-", (age_group*5-1))]
  
  ref_breaks <- c(1, 5, 10, 20, 2000)
  categ <- c("1-5","5-10","10-20","20+")
  Total <- Total %>% 
    mutate(category=cut(clusters, breaks=ref_breaks,include.lowest = T, 
                        labels=categ))
  
  p <- ggplot(Total, aes(x=Total$long, y=Total$lat))
  p <- p + geom_polygon(data=Total, aes(group = group, fill = Total$category),
                        color="black") +
    scale_fill_manual(values = brewer.pal(n = 4, name = 'Purples'), 
                      na.value = "grey", breaks = categ)
  
  P1 <- p + theme_bw()  + labs(fill = "", x="", y="", title = "A")
  P1 <- P1 + scale_y_continuous(breaks=c()) + scale_x_continuous(breaks=c()) + 
    theme(panel.border =  element_blank(), legend.text=element_text(size=20), 
          legend.title = element_text(size = 18), title = element_text(size = 20))
  
  P2 <- ggplot(dt_cases, 
               aes(x = factor(dt_cases$age_group_char, 
                              levels = paste0(((1:max(dt_cases$age_group))-1)*5, "-",
                                              ((1:max(dt_cases$age_group))*5-1)))))+
    geom_histogram(colour = NA, stat = "count")
  P2 <- P2 + theme_classic() + theme(axis.title = element_text(size = 18),
                                     axis.text = element_text(size = 15), 
                                     title = element_text(size = 20))  + 
    labs(x="Age group", y="Number of cases", title = "B")
  
  grid.arrange(P1, P2, nrow = 2)
  
}

supp_post <- function(list_out, burnin, sample_every){
  par(mfrow = c(length(list_out),1), oma = c(2,1,0,1), cex.lab = 2, 
      cex.axis = 1.5, mar = c(1,5,1,-1)+ 1.1, las=1, bty="l", tcl=-1, lwd=2,
      mgp = c(3, 1, 0))
  for(i in 1:length(list_out)){
    out_X <- list_out[[i]]
    out_X_burnin <- as.numeric(out_X[(burnin/sample_every):dim(out_X)[1],
                                     "post"])
    plot(out_X_burnin, type = "l", ylab = "")
    mtext(text = LETTERS[i], adj = 0, cex = 1.5, line = 0.5)
  }
  title(outer = T, ylab = "Posterior", line = -1, cex.lab = 2)
  title(outer = T, xlab = "Iteration", line = 0.5, cex.lab = 2)
  
}
