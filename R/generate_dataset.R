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
  
  genotype_k <- rep(genotype_index, sec)
  
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
         r0_state, pop_county, t_min, t_max){
  
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
    date_i <- sample(x = t_min:t_max, size = 1)
    date_i <- as.Date(date_i, origin = "1970-01-01")
    
    county_i <- sample(x =names(pop_county), size = 1, 
                     prob = pop_county/sum(pop_county))
    
    age_group_i <- ceiling(rexp(n = 1, rate = 0.2))
    if(age_group_i>dim(polymod_prop)[2]) age_group_i <- dim(polymod_prop)[2]
    
    genotype_i <- sample(c("B3", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "G3", "H1", "H2"),
                         size = 1, 
                         prob = c(0.36, 0.05, 0.10, 0.05, 0.01, 0.01, 0.25, 0.10, 
                                  0.01, 0.05, 0.01))
    
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
  
  gen_reported <- runif(nrow(dt_cases), 0, 1)
  dt_cases[gen_reported < 0.6, Genotype := "Not attributed"]
  
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