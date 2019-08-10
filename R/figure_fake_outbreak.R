#### Figure 2 ####

png(filename = paste0(data.out, "Figure_2_fake.png"), width = 1000, 
    height = 750)
#Create figure window and layout
plot.new()
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))

dt_heatmap <- rbind(dt_heatmap_001, dt_heatmap_005, dt_heatmap_import)

pushViewport(viewport(layout.pos.row = 1))
par(fig = gridFIG(),family = "sans",
    omi = c(0,0,0,0),
    cex.lab = 1.4,
    cex.axis = 1.2, 
    mar = c(3,3.5,0,1)+ 1.1,
    las=1, bty="l", tcl=-1, lwd=2,
    mgp = c(3, 1, 0), new = T)

b <- barplot(rbind(med_size_cluster_barplot_001[-1],
                   med_size_cluster_barplot_005[-1],
                   med_size_cluster_barplot_import[-1],
                   size_cluster_ref_barplot[-1]), 
             col = grey.colors(6)[c(2, 3, 4, 6)],
             beside = T, 
             ylim = c(0, max(size_cluster_ref_barplot[-1],
                             up_size_cluster_barplot_001[-1],
                             up_size_cluster_barplot_005[-1],
                             up_size_cluster_barplot_import[-1])), 
             xlab = "Cluster size", ylab = "Proportion", 
             border = NA)

arrows(x0 = b[1, ], y0 = low_size_cluster_barplot_001[-1], 
       x1 = , y1 = up_size_cluster_barplot_001[-1], angle = 90, 
       code = 3,length = 0.1, lwd =  0.5)
arrows(x0 = b[2, ], y0 = low_size_cluster_barplot_005[-1], 
       x1 = , y1 = up_size_cluster_barplot_005[-1], angle = 90, 
       code = 3,length = 0.1, lwd =  0.5)
arrows(x0 = b[3, ], y0 = low_size_cluster_barplot_import[-1], 
       x1 = , y1 = up_size_cluster_barplot_import[-1], angle = 90, 
       code = 3,length = 0.1, lwd =  0.5)

legend("right",
       fill = c(grey.colors(6)[c(2, 3, 4, 6)], NA), 
       cex = 1.4,
       border = NA, 
       pch = c(NA, NA, NA, NA, "_"), pt.cex = 3,
       legend = c("No import, threshold 1%", 
                  "No import, threshold 5%", 
                  "With import",
                  "Correct inference"), bty = "n")

par(fig = c(0.37, 0.78, 0.72, 1),family = "sans",
    omi = c(0,0,0,0),
    cex.lab = 1.4,
    cex.axis = 1.2, 
    mar = c(3,3.5,0,1)+ 1.1,
    las=1, bty="l", tcl=-1, lwd=2,
    mgp = c(3, 1, 0), new = T)

nb_clust_ref <- length(unique(dt_cases$cluster))
names(nb_clust_ref) <- "Imports"
b <- barplot(rbind(c(med_001[1], 
                     Singletons = as.numeric(med_size_cluster_barplot_001[1])) 
                   , c(med_005[1], 
                       Singletons = as.numeric(med_size_cluster_barplot_005[1])) 
                   , c(med_import[1], 
                       Singletons = as.numeric(med_size_cluster_barplot_import[1])) 
                   , c(nb_clust_ref, 
                       Singletons = as.numeric(size_cluster_ref_singletons[1]))), 
             col = grey.colors(6)[c(2, 3, 4, 6)], beside = T, 
             ylim = c(0, max(c(up_001[1], up_005[1], up_import[1],
                               nb_clust_ref))), ylab = "Number of cases", 
             border = NA)

arrows(x0 = b[1, ], y0 = c(low_001[1], low_size_cluster_barplot_001[1]), 
       x1 = , y1 = c(up_001[1], up_size_cluster_barplot_001[1]), angle = 90, 
       code = 3,length = 0.1, lwd =  0.5)
arrows(x0 = b[2, ], y0 = c(low_005[1], low_size_cluster_barplot_005[1]), 
       x1 = , y1 = c(up_005[1], up_size_cluster_barplot_005[1]), angle = 90, 
       code = 3,length = 0.1, lwd =  0.5)
arrows(x0 = b[3, ], y0 = c(low_import[1], low_size_cluster_barplot_import[1]), 
       x1 = , y1 = c(up_import[1], up_size_cluster_barplot_import[1]), angle = 90, 
       code = 3,length = 0.1, lwd =  0.5)

sing_vect <- c(rbind(c(med_001[2], med_prop_singletons_001[1]),
                     c(med_005[2], med_prop_singletons_005[1]),
                     c(med_import[2], med_prop_singletons_import[1]),
                     c(size_cluster_ref_barplot %>% sum, size_cluster_ref_singletons[1])))

points(sing_vect[sing_vect>0] ~ c(b)[sing_vect > 0], pch = "_", cex = 3)

popViewport()

pushViewport(viewport(layout.pos.row = 2))

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

dev.off()

#### Figure 4 ####
