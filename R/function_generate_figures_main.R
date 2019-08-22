generate_figure_3 <- function(dt_cases, ref_breaks, categ){
  data(state, package = "datasets")
  dt_map_cases <- data.table(abb = c(state.abb, "COL"))
  dt_map_cases[, region := tolower(c(state.name, "district of columbia"))]
  setkey(dt_map_cases, abb)
  dt_map_cases[, cases := 0]
  dt_map_cases[names(table(dt_cases$State)),
               cases := as.numeric(table(dt_cases$State))]
  dt_map_cases <- dt_map_cases[region != "alaska",]
  dt_map_cases <- dt_map_cases[region != "hawaii",]
  
  all_states <- map_data("state")
  Total <- merge(all_states, dt_map_cases, by="region")
  Total <- Total %>% mutate(category=cut(cases, breaks=ref_breaks,
                                         include.lowest = T, labels=categ))
  
  p <- ggplot(Total, aes(x=Total$long, y=Total$lat))
  p <- p + geom_polygon(data=Total, aes(group = group, fill = Total$category),
                        color="black") +
    scale_fill_manual(values = brewer.pal(n = length(ref_breaks), name = 'Purples'),
                      na.value = "grey", breaks = categ)

  P1 <- p + theme_bw()  + labs(fill = "Number of \ncases", fill = "", x="", y="",
                               title = "A")
  P1 <- P1 + scale_y_continuous(breaks=c()) + scale_x_continuous(breaks=c()) + 
    theme(panel.border =  element_blank(), legend.text=element_text(size=20), 
          legend.title = element_text(size = 18), title = element_text(size = 20))
  
  P2 <- ggplot(dt_cases, aes(x = year(dt_cases$Date))) +
    geom_histogram(colour = NA, binwidth = 1)
  P2 <- P2 + theme_classic() + theme(axis.title = element_text(size = 18),
                                     axis.text = element_text(size = 15), 
                                     title = element_text(size = 20))  + 
    labs(x="Year", y="Number of cases", title = "B")
  
  grid.arrange(P1, P2, nrow = 2)
  
}

#' Title: Generate Figure 3 / 4
#'
#'
#' @param fig_hist_list list of output list for prepare_for_figure(), for every
#' scenario to be displayed on the histogram
#' @param list_fig_heatmap list of output list for prepare_for_figure(), for every
#' scenario to be displayed on the heatmaps
#'
#' @export
#'
#' @examples
generate_figure_4_5 <- function(fig_hist_list, list_fig_heatmap, ref = T){
  ## Median values of cluster size distributions
  med_size_cluster_list <- lapply(fig_hist_list,
                                  function(X)
                                    if(!is.null(X$med_size_cluster_barplot))
                                      return(as.data.frame(t(X$med_size_cluster_barplot))))
  med_size_cluster <- bind_rows(med_size_cluster_list) %>% as.matrix
  ## 97.5% CI
  up_size_cluster_list <- lapply(fig_hist_list,
                                  function(X)
                                    if(!is.null(X$up_size_cluster_barplot))
                                    return(as.data.frame(t(X$up_size_cluster_barplot))))
  up_size_cluster <- bind_rows(up_size_cluster_list) %>% as.matrix
  
  ## 2.5% CI
  low_size_cluster_list <- lapply(fig_hist_list,
                                  function(X)
                                    if(!is.null(X$low_size_cluster_barplot))
                                      return(as.data.frame(t(X$low_size_cluster_barplot))))
  low_size_cluster <- bind_rows(low_size_cluster_list) %>% as.matrix
  
  ## Number of singletons
  singletons_list <- lapply(fig_hist_list,
                                  function(X)
                                    if(!is.null(X$med_prop_singletons))
                                      return(as.data.frame(t(X$med_prop_singletons))))
  singletons <- bind_rows(singletons_list) %>% as.matrix
  
  ## Median number of imports, and number of imports correctly inferred
  med_imports_list <- lapply(fig_hist_list,
                            function(X)
                              if(!is.null(X$med))
                                return(as.data.frame(t(X$med))))
  med_imports <- bind_rows(med_imports_list) %>% as.matrix

  ## 97.5% CI
  up_imports_list <- lapply(fig_hist_list,
                             function(X)
                               if(!is.null(X$up))
                                 return(as.data.frame(t(X$up))))
  up_imports <- bind_rows(up_imports_list) %>% as.matrix
  
  ## 2.5% CI
  low_imports_list <- lapply(fig_hist_list,
                             function(X)
                               if(!is.null(X$low))
                                 return(as.data.frame(t(X$low))))
  low_imports <- bind_rows(low_imports_list) %>% as.matrix

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
  if(ref == F){
    arrows(x0 = b, y0 = low_size_cluster[, -1], x1 = , 
           y1 = up_size_cluster[, -1], angle = 90, code = 3,
           length = 0.1, lwd =  0.5)
  } else
    arrows(x0 = b[-nrow(b),], y0 = low_size_cluster[, -1], x1 = , 
           y1 = up_size_cluster[, -1], angle = 90, code = 3,
           length = 0.1, lwd =  0.5)

  
  ## Number of singletons and imports
  par(fig = c(0.27, 0.72, 0.72, 1),family = "sans",
      omi = c(0,0,0,0), cex.lab = 1.4, cex.axis = 1.2, 
      mar = c(3,3.5,0,1)+ 1.1, las=1, bty="l", tcl=-1, lwd=2,
      mgp = c(3, 1, 0), new = T)
  b <- barplot(cbind(imports = med_imports[,1], singletons = med_size_cluster[, 1]),
               col = grey.colors(nrow(med_size_cluster)),
               beside = T,
               ylim = c(0, max(cbind(up_imports[,1], up_size_cluster[, 1]))), 
               ylab = "Number of cases", 
               border = NA)
  # CIs
  if(ref == F){
    arrows(x0 = b, y0 = cbind(low_imports[,1], low_size_cluster[, 1]), 
           x1 = , y1 = cbind(up_imports[,1], up_size_cluster[, 1]), angle = 90, 
           code = 3,length = 0.1, lwd =  0.5)
  } else 
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
  
  ## List of all heatmaps in data
  list_heatmap <- lapply(list_fig_heatmap, function(X) return(X$dt_heatmap))
  ## dt_heatmap merges the data table contained in list_heatmap
  dt_heatmap <- bind_rows(list_heatmap, .id = "id")
  dt_heatmap[prop<0.001, prop := 0]
  # Generate the heatmaps
  p <- ggplot(dt_heatmap, aes(sensitivity, precision)) + 
    geom_tile(aes(fill = prop)) + 
    facet_grid(.~id) + 
    scale_fill_gradient2(na.value = "lightgrey",
                         low = "white", mid = "lightblue", midpoint = log(0.1),
                         high = "darkblue",
                         limits = c(0.005,1),
                         trans = "log", breaks = c(0,0.01, 0.1,0.5),
                         name = "Proportion 
of cases")
  p <- p + theme_classic() + 
    labs(x = "Sensitivity",
         y = "Precision") +
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
generate_figure_6 <- function(list_factor_import, ref_breaks, categ){
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