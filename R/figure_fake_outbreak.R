setwd(dir = "/PhD/Study/Measles_clustering_susc/")

library(data.table)
library(measlesoutbreaker)
library(geosphere)
library(xlsx)
library(dplyr)
library(Rcpp)
library(grid)
library(ggplot2)
library(gridBase)
library(gridExtra)


data("fake_outbreak")
source("src/prepare_out_figures.R")

dt_cases <- fake_outbreak[["cases"]]
dist_mat <- fake_outbreak[["distance"]]
pop_vect <- fake_outbreak[["population"]]
age_contact <- fake_outbreak[["age_contact"]]
burnin <- 10000
sample_every <- 50
dt_cases <- dt_cases[order(Date), ]
max_clust <- 500

dt_cases[, size_cluster := table(cluster)[cluster]]

thresh_barplot <- c(1,2,3,5,10,15,20,30,40,50)

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

diff <- 0.25
ref <- seq(diff/2, 1 - diff/2, diff)
missed <- seq(diff/2, 1 - diff/2, diff)


#### Reference clusters ####

ref_size_clusters <- dt_cases[!duplicated(cluster) & cluster != ".", size_cluster]
ref_size_clusters <- c(rep(1, dim(dt_cases[cluster == ".",])[1]), ref_size_clusters)
ref_size_clusters <- factor(ref_size_clusters, levels = 1:max_clust)
ref_size_clusters <- table(ref_size_clusters)
size_cluster_ref_barplot <- t(rowsum(x = as.matrix(ref_size_clusters), 
                                     group = groups_barplot))
colnames(size_cluster_ref_barplot) <- names_barplot

size_cluster_ref_singletons <- numeric(length(size_cluster_ref_barplot))
names(size_cluster_ref_singletons) <- names(size_cluster_ref_barplot)
size_cluster_ref_singletons[1] <- size_cluster_ref_barplot[1]

#### threshold 001 ####

out <- readRDS(file = "fake_outbreak_runs/no_import_thresh001.rds")
fig_001 <- prepare_for_figures(out = out, 
                               dt_cases = dt_cases)

low_size_cluster_barplot_001 <- fig_001[["low_size_cluster_barplot"]]
med_size_cluster_barplot_001 <- fig_001[["med_size_cluster_barplot"]]
up_size_cluster_barplot_001 <- fig_001[["up_size_cluster_barplot"]]
med_prop_singletons_001 <- fig_001[["med_prop_singletons"]]
low_001 <- fig_001[["low"]]
med_001 <- fig_001[["med"]]
up_001 <- fig_001[["up"]]
dt_heatmap_001 <- fig_001[["dt_heatmap"]]
dt_heatmap_001[, type := "001"]

#### threshold 005 ####

out <- readRDS(file = "fake_outbreak_runs/no_import_thresh005.rds")

fig_005 <- prepare_for_figures(out = out, 
                               dt_cases = dt_cases)
low_size_cluster_barplot_005 <- fig_005[["low_size_cluster_barplot"]]
med_size_cluster_barplot_005 <- fig_005[["med_size_cluster_barplot"]]
up_size_cluster_barplot_005 <- fig_005[["up_size_cluster_barplot"]]
med_prop_singletons_005 <- fig_005[["med_prop_singletons"]]
low_005 <- fig_005[["low"]]
med_005 <- fig_005[["med"]]
up_005 <- fig_005[["up"]]
dt_heatmap_005 <- fig_005[["dt_heatmap"]]
dt_heatmap_005[, type := "005"]

#### with imports ####

out <- readRDS(file = "fake_outbreak_runs/with_import.rds")

fig_import <- prepare_for_figures(out = out, 
                                  dt_cases = dt_cases)
low_size_cluster_barplot_import <- fig_import[["low_size_cluster_barplot"]]
med_size_cluster_barplot_import <- fig_import[["med_size_cluster_barplot"]]
up_size_cluster_barplot_import <- fig_import[["up_size_cluster_barplot"]]
med_prop_singletons_import <- fig_import[["med_prop_singletons"]]
low_import <- fig_import[["low"]]
med_import <- fig_import[["med"]]
up_import <- fig_import[["up"]]
dt_heatmap_import <- fig_import[["dt_heatmap"]]
dt_heatmap_import[, type := "import"]


#### Figure 2 ####

data.out <- "H:/PhD/Study/Measles_clustering_susc/fake_outbreak_runs/"

png(filename = paste0(data.out, "Figure_2_fake.png"), width = 1000, 
    height = 750)
#Create figure window and layout
plot.new()
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 1)))

dt_heatmap <- rbind(dt_heatmap_001, dt_heatmap_005, dt_heatmap_import)
# dt_heatmap <- dt_heatmap_005

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
