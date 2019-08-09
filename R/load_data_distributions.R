source("R/library_importation.R")

#### Serial interval Measles ####

# Serial interval: 
# Normal distribution, mean = 11.7, sd = 2.0 (source : Vink et al 2014)
# w_rd <- rnorm(1e6, mean = 11.7, sd = 2.0)
# hist(w_rd, breaks = seq(min(w_rd), max(w_rd)+1, 1))
# w <- hist(w_rd, freq = F, breaks=seq(trunc(min(w_rd)),
#                                      max(ceiling(w_rd))))$de

w <- dnorm(x = 1:300, mean = 11.7, sd = 2.0, log = T)

#### Incubation period ####

#From "The correlation between infectivity and incubation period of measles, 
# estimated from households with two cases"
# Don Klinkenberg, Hiroshi Nishiura
# mean = shape * scale = 11.5
# variance = shape * scale**2 = 2.2**2
# scale = variance / mean = 0.4286
# shape = mean / scale = 26.834
# f_rd <-rgamma(n = 1e6, scale = 0.4286, shape = 26.834)
# hist(f_rd, freq = F, breaks=seq(trunc(min(f_rd)),
#                                      max(ceiling(f_rd))))
f <- dgamma(x = 1:300, scale = 0.4286, shape = 26.834, log = T)



#### Contact matrices ####

# From Polymod
polymod <- as.data.table(read.csv(file = "data/social_contact_UK.csv",
                                  header = T, , colClasses = "numeric"))
colnames(polymod) <- substr(x = colnames(polymod), start = 2, stop = nchar(colnames(polymod)))
#row <- age of contact
#column <- age of participant
polymod_prop <- apply(polymod, 2, function(X) 
  return(X/sum(X)))


#### Distance Counties pop ####

dt_loca_pop_center <- as.data.table(read.csv(
  file = "data/pop_center.csv",
  stringsAsFactors = FALSE, colClasses = c("character", "character",
                                           "character","character",
                                           "numeric", "numeric",
                                           "numeric", "character")))
dt_state_county <- dt_loca_pop_center[,.(ID_COUNTY, STNAME)]
setkey(dt_state_county, ID_COUNTY)
dt_distance_pop_center <- as.data.table(matrix(0, nrow = nrow(dt_loca_pop_center)**2,
                                               ncol = 9))
colnames(dt_distance_pop_center) <- c("county1", "county2", "distance_km",
                                      "pop_county1", "long1", "lat1",
                                      "pop_county2", "long2", "lat2")
dt_distance_pop_center[, county1 := as.character(county1)]
dt_distance_pop_center[, county1 := rep(dt_loca_pop_center$ID_COUNTY, 
                                        nrow(dt_loca_pop_center))]
setkey(dt_loca_pop_center, ID_COUNTY)
dt_distance_pop_center[, long1 := dt_loca_pop_center[dt_distance_pop_center$county1,
                                                     LONGITUDE]]
dt_distance_pop_center[, lat1 := dt_loca_pop_center[dt_distance_pop_center$county1,
                                                    LATITUDE]]
dt_distance_pop_center[, county2 := rep(dt_loca_pop_center$ID_COUNTY,
                                        each = nrow(dt_loca_pop_center))]
dt_distance_pop_center[, long2 := dt_loca_pop_center[dt_distance_pop_center$county2,
                                                     LONGITUDE]]
dt_distance_pop_center[, lat2 := dt_loca_pop_center[dt_distance_pop_center$county2,
                                                    LATITUDE]]
long1 <- dt_distance_pop_center$long1
lat1 <- dt_distance_pop_center$lat1
long2 <- dt_distance_pop_center$long2
lat2 <- dt_distance_pop_center$lat2
mat1 <- matrix(c(long1,
                 lat1), ncol = 2)
mat2 <- matrix(c(long2,
                 lat2), ncol = 2)
rm(list = c("lat1", "lat2", "long1", "long2"))
dist <- numeric(dim(dt_distance_pop_center)[1])
gc()
dist <- distGeo(mat1, mat2)/1000
dt_distance_pop_center[, distance_km := dist]

dt_distance_pop_center[, pop_county1 := 
                         dt_loca_pop_center[dt_distance_pop_center$county1,
                                            POPULATION]]
dt_distance_pop_center[, pop_county2 := 
                         dt_loca_pop_center[dt_distance_pop_center$county2,
                                            POPULATION]]
dt_distance <- dt_distance_pop_center

rm(list = c("dist", "mat1", "mat2", "dt_distance_pop_center", "dt_loca_pop_center"))
gc()

dt_population <- dt_distance[!duplicated(county1), .(county1, pop_county1)]
colnames(dt_population) <- c("county", "population")


