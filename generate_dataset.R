#### Libraries ####

library(magrittr)
library(dplyr)
library(expm)
library(devtools)
library(Rcpp)
library(ggplot2)
library(data.table)
library(geosphere)

#### Serial interval Measles ####

# Serial interval: 
# Normal distribution, mean = 11.7, sd = 2.0 (source : Vink et al 2014)
# w_rd <- rnorm(1e6, mean = 11.7, sd = 2.0)
# hist(w_rd, breaks = seq(min(w_rd), max(w_rd)+1, 1))
# w <- hist(w_rd, freq = F, breaks=seq(trunc(min(w_rd)),
#                                      max(ceiling(w_rd))))$de

w <- dnorm(x = 1:300, mean = 11.7, sd = 2.0, log = T)

#### Incubation period ####

#From The correlation between infectivity and incubation period of measles, estimated from households with two cases
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


#### Function to create new cases ####

new_case <- function(sec, date_index, county_index, age_group_index, genotype_index, 
                     count, cluster, generation, dt_distance, 
                     polymod_prop, w){
  county_k <- sample(dt_distance[county2 == county_index, county1], sec, 
                     prob = dt_distance[county2 == county_index, probs], 
                     replace = T)
  
  age_group_k <- sample(1:15, size = sec, prob = polymod_prop[, age_group_index],
                        replace = T)
  
  date_k <- date_index+sample(1:length(w), size = sec, prob = exp(w),
                              replace = T)
  
  genotype_k <- sample(c("Not attributed", genotype_index), size = sec, 
                       prob = c(0.60, 0.40), replace = T)
  case_k <- data.table(as.character((count):(count + sec - 1)), 
                       dt_state_county[county_k, STNAME],
                       date_k, 
                       genotype_k, 
                       county_k, 
                       age_group_k, 
                       rep(FALSE, sec), 
                       cluster, 
                       generation + 1)
  colnames(case_k) <- c("ID", "State", "Date", "Genotype", "county", "age_group", 
                        "import", "cluster", "generation")
  
  return(case_k)
}

#### Generate dataset ####


generate_dataset <- function(a, b, gamma, dt_distance, polymod_prop, w, nb_cases, 
         r0_state, pop_county){
  
  dt_distance[, nb_commut := pop_county1**a*exp(-distance_km*b)]
  dt_distance[distance_km > gamma, nb_commut := 0]
  sum_commut <- dt_distance[, lapply(.SD,sum), by = county2,
                            .SDcols = "nb_commut"]
  setkey(sum_commut, county2)
  dt_distance[, probs := nb_commut / sum_commut[county2,nb_commut]]
  
  dt_cases <- data.table(NULL)
  count <- 1
  i <- 1
  
  while(dim(dt_cases)[1] < nb_cases){
    date_i <- sample(x = as.Date("2010-01-01"):as.Date("2016-01-01"), size = 1)
    date_i <- as.Date(date_i, origin = "1970-01-01")
    
    county_i <- sample(x =names(pop_county), size = 1, 
                     prob = pop_county/sum(pop_county))
    
    age_group_i <- ceiling(rexp(n = 1, rate = 0.2))
    if(age_group_i>dim(polymod_prop)[2]) age_group_i <- dim(polymod_prop)[2]
    
    genotype_i <- sample(c("Not attributed", "D4", "B1", "B4", "D8"), size = 1, 
                         prob = c(0.60, 0.20, 0.10, 0.05, 0.05))
    
    clust_i <- i
    gen <- 1
    case_i <- data.table(as.character(count), dt_state_county[county_i, STNAME],
                         date_i, genotype_i, county_i, age_group_i, TRUE, clust_i,
                         gen)
    colnames(case_i) <- c("ID", "State", "Date", "Genotype", "county", 
                          "age_group", "import", "cluster", "generation")
    dt_cases <- rbind(dt_cases, case_i)
    count <- dim(dt_cases)[1]
    while(count<= dim(dt_cases)[1]
          ){
      date_i <- dt_cases[count, Date]
      county_i <- dt_cases[count, county]
      age_group_i <- dt_cases[count, age_group]
      genotype_i <- dt_cases[count, Genotype]
      clust_i <- dt_cases[count, cluster]
      r0_i <- r0_state[dt_state_county[county_i, STNAME]]
      gen_i <- dt_cases[count, generation]
      new_cases_i <- rpois(1, r0_i)
      count <- count + 1
      if(new_cases_i > 0){
        case_subs <- new_case(sec = new_cases_i, date_index = date_i, 
                              county_index = county_i, age_group_index = age_group_i, 
                              genotype_index = genotype_i, count = count, 
                              cluster = clust_i, generation = gen_i, 
                              dt_distance = dt_distance,
                              polymod_prop = polymod_prop, w = w)
        dt_cases <- rbind(dt_cases, case_subs)
      }
    }
    i <- i + 1 
  }
  
  dt_distance <- dt_distance[,.(county1, county2, distance_km)]
  
  setkey(dt_distance, county1)
  distance_matrix <- sapply(unique(dt_distance$county1), 
                            function(X){
    return(dt_distance[county2 == X][, distance_km])
  })
  rownames(distance_matrix) <- colnames(distance_matrix)
  
  vect_pop <- dt_population$population
  names(vect_pop) <- dt_population$county
  vect_pop <- vect_pop[rownames(distance_matrix)]
  
  fake_outbreak <- list(cases = dt_cases, distance = distance_matrix,
                        population = vect_pop, age_contact = polymod)
  return(fake_outbreak)
}