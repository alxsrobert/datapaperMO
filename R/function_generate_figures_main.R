#' Title: Generate Figure 3 / 4
#'
#' @param med_size_cluster List containing the Median values of cluster size 
#' distributions in each run or data included in the histogram.
#' @param up_size_cluster List containing the 97.5% quantile of cluster size 
#' distributions in each run or data included in the histogram.
#' @param low_size_cluster List containing the 2.5% quantile of cluster size 
#' distributions in each run or data included in the histogram.
#' @param singletons List containing the 97.5% quantile of cluster size 
#' distributions in each run or data included in the histogram.
#' @param med_imports List containing the median number of imports, and the
#' number of imports correctly inferred in each run or data included in the
#' histogram.
#' @param up_imports List containing the 97.5% quantile of the number of imports, 
#' and the number of imports correctly inferred in each run or data included in 
#' the histogram.
#' @param low_imports List containing the 97.5% quantile of the number of imports, 
#' and the number of imports correctly inferred in each run or data included in the
#' histogram.
#' @param list_heatmap List containing the heatmap data tables for each run.
#'
#' @export
#'
#' @examples
generate_figure_3_4 <- function(med_size_cluster, up_size_cluster, 
                                low_size_cluster, singletons, med_imports,
                                up_imports, low_imports,list_heatmap){
  ## Plot is mixing base plot and ggplot so need to use Viewport
  plot.new()
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))
  
  ## First: Histograms
  pushViewport(viewport(layout.pos.row = 1))
  par(fig = gridFIG(),family = "sans",
      omi = c(0,0,0,0),
      cex.lab = 1.4,
      cex.axis = 1.2, 
      mar = c(3,3.5,0,1)+ 1.1,
      las=1, bty="l", tcl=-1, lwd=2,
      mgp = c(3, 1, 0), new = T)
  # Cluster size distribution in each run, except singletons
  b <- barplot(med_size_cluster[, -1], 
               col = grey.colors(nrow(med_size_cluster)),
               beside = T, 
               ylim = c(0, max(up_size_cluster[, -1])), 
               xlab = "Cluster size", ylab = "Proportion", 
               border = NA)
  # CIs
  arrows(x0 = b[-nrow(b),], y0 = low_size_cluster[, -1], 
         x1 = , y1 = up_size_cluster[, -1], angle = 90, 
         code = 3,length = 0.1, lwd =  0.5)

  
  ## Number of singletons and imports
  par(fig = c(0.37, 0.78, 0.72, 1),family = "sans",
      omi = c(0,0,0,0),
      cex.lab = 1.4,
      cex.axis = 1.2, 
      mar = c(3,3.5,0,1)+ 1.1,
      las=1, bty="l", tcl=-1, lwd=2,
      mgp = c(3, 1, 0), new = T)
  b <- barplot(cbind(imports = med_imports[,1], singletons = med_size_cluster[, 1]),
               col = grey.colors(nrow(med_size_cluster)),
               beside = T,
               ylim = c(0, max(cbind(up_imports[,1], up_size_cluster[, 1]))), 
               ylab = "Number of cases", 
               border = NA)
  # CIs
  arrows(x0 = b[-nrow(b), ], y0 = cbind(low_imports[,1], low_size_cluster[, 1]), 
         x1 = , y1 = cbind(up_imports[,1], up_size_cluster[, 1]), angle = 90, 
         code = 3,length = 0.1, lwd =  0.5)
  # Number of imports / singletons correctly inferred (i.e. who are 
  # imports / singletons in the data)
  sing_vect <- c(med_imports[,2], singletons[,1])

  points(sing_vect[sing_vect>0] ~ c(b)[sing_vect > 0], pch = "_", cex = 3)
  
  popViewport()
  
  ## Heatmap
  pushViewport(viewport(layout.pos.row = 2))
  
  ## dt_heatmap merges the data table contained in list_heatmap
  dt_heatmap <- bind_rows(list_heatmap)
  # Generate the heatmaps
  p <- ggplot(dt_heatmap, aes(ref, missed)) + 
    geom_tile(aes(fill = prop)) + 
    facet_grid(.~type) + 
    scale_fill_gradient2(na.value = "lightgrey",
                         low = "white", mid = "lightblue", midpoint = log(0.1),
                         high = "darkblue",
                         limits = c(0.005,1),
                         trans = "log", breaks = c(0,0.01, 0.1,0.5),
                         name = "Proportion 
of cases")
  p <- p + theme_classic() + labs(x = "Proportion of epi cluster in inferred cluster",
                                  y = "Proportion of inferred cluster not in epi cluster
(False positives)") +
    theme(axis.ticks = element_blank(), 
          axis.title = element_text(size = rel(1.4)),
          legend.text = element_text(size = rel(1.2)), 
          legend.title = element_text(size = rel(1.4)),
          axis.text = element_text(size = rel(1.2)),
          strip.text.x = element_blank()
          , text = element_text(family = "")
    ) + xlim(0,1) + ylim(0,1)
  print(p, newpage = FALSE)
  popViewport()
}

#' Title: Generate Figure 5
#'
#' @param list_factor_import 
#' @param ref_breaks 
#' @param categ 
#'
#' @export
#'
#' @examples
generate_figure_5 <- function(list_factor_import, ref_breaks, categ){
  ## Load data with all abbreviations and state names
  data(state, package = "datasets")
  dt_map_cases <- data.table(abb = c(state.abb, "COL"))
  dt_map_cases[, region := tolower(c(state.name, "district of columbia"))]
  setkey(dt_map_cases, abb)
  dt_map_cases[, cases := 0]
  
  ## Merge the state factors obtained in each run in a list, then in a data table
  list_map_cases <- lapply(list_factor_import, function(X){
    dt_map_cases[colnames(X), fact := apply(X, 2, median)]
    dt_map_cases <- dt_map_cases %>% 
      mutate(category=cut(fact, breaks=ref_breaks,include.lowest = T,
                          labels=categ))
    return(dt_map_cases)
  })
  dt_map <- bind_rows(list_map_cases, .id = "type")

  all_states <- map_data("state")
  Total <- merge(all_states, dt_map, by="region")
  p <- ggplot(Total, aes(x=Total$long, y=Total$lat)) +   
    facet_grid(factor(type)~.)
  p <- p + geom_polygon(data=Total, aes(group = group
                                        , fill = Total$category),
                        color = "black") +
    scale_fill_manual(values = brewer.pal(n = 5, name = 'Purples'),
                      na.value = "grey"
                      , breaks = categ
    )
  p <- p + theme_bw()  + 
    labs(fill = "", x="", y="")
  p <- p + scale_y_continuous(breaks=c()) + scale_x_continuous(breaks=c()) + 
    theme(panel.border =  element_blank(), legend.text=element_text(size=25), 
          strip.text.y = element_blank(),
          legend.title = element_text(size = 18), title = element_text(size = 20))
  
  print(p)

}