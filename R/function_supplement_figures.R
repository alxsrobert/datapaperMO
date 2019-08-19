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

estim_params <- function(out){
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

supp_fig_stratified_state <- function(out){
  
}
supp_fig_sec_overall <- function(list_out){
  
}
supp_fig_sec_map <- function(list_out){
  
}

supp_fig_distance_transmission <- function(med_dist, low_dist, up_dist){
  par(mfrow = c(1,1), family = "sans", omi = c(0,0,0,0), cex.lab = 1.4,
      cex.axis = 1.2,  mar = c(3,3.5,0,1)+ 1.1, las=1, bty="l", tcl=-1, lwd=2,
      mgp = c(3, 1, 0))
  
  b <- barplot(rbind(med_dist), 
               col =grey.colors(nrow(med_size_cluster)),
               beside = T, 
               ylim = c(0, 1),
               xlab = "Distance (km)", ylab = "Proportion", border = NA)
  arrows(x0 = b[1, ], y0 = low_dist, x1 = , y1 = up_dist, angle = 90, code = 3,
         length = 0.1)

}