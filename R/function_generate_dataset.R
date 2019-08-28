#' Title: Function to create new cases
#'
#' @param sec Number of new cases to generate.
#' @param date_index Infecion date of the index case.
#' @param county_index County of the index case.
#' @param age_group_index Age group of the index case.
#' @param genotype_index Genotype of the index case.
#' @param count ID of the index case.
#' @param cluster Cluster of the index case.
#' @param generation Generation number of the index case.
#' @param dt_distance Data table of the distance between counties.
#' @param polymod_prop Matrix of proportion of contacts between age groups.
#' @param w Vector of the serial interval of the disease.
#' @param dt_state_county Data table of each county and the state they belong to.
#' @param prop_gen: proportion of genotype set to NA after generating the data
#'
#' @return
#' 
#' A data frame containing the ID, State, Date, Genotype, County, Age group, 
#' import status, cluster and generation number of the generated cases
#' @export
#'
#' @examples
new_case <- function(sec, date_index, county_index, age_group_index, genotype_index, 
                     count, cluster, generation, dt_distance, 
                     polymod_prop, w, dt_state_county){
  # Draw new cases' county
  county_k <- sample(dt_distance[county2 == county_index, county1], sec, 
                     prob = dt_distance[county2 == county_index, probs], 
                     replace = T)
  # Draw new cases' age
  age_group_k <- sample(1:15, size = sec, prob = polymod_prop[, age_group_index],
                        replace = T)
  # Draw new cases' infection date
  date_k <- date_index+sample(1:length(w), size = sec, prob = exp(w),
                              replace = T)
  # new cases' genotype is the same as index's
  genotype_k <- rep(genotype_index, sec)
  # Merge generated features into a data table
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

#' Title: Generate dataset
#'
#' @param a double first spatial parameter (population).
#' @param b double second spatial parameter (distance).
#' @param gamma double distance threshold.
#' @param dt_distance data table of the distance between counties.
#' @param polymod_prop Matrix of proportion of contacts between age groups.
#' @param w Vector of the serial interval of the disease.
#' @param nb_cases integer maximum number of cases to be generated.
#' @param r0_state vector r0 in each state.
#' @param pop_county
#' @param t_min Date minimum date of infection of ancestors.
#' @param t_max Date maximum date of infection of ancestors.
#' @param dt_state_county Data table of each county and the state they belong to.
#' @param prop_gen: proportion of genotype set to NA after generating the data
#'
#' @return
#' A list containing the linelist of cases, and the distance matrix, the number of 
#' inhabitants in each county, and the age matrix used to generate it.
#' @export
#'
#' @examples
generate_dataset <- function(a, b, gamma, dt_distance, 
                             polymod_prop, w, nb_cases, r0_state, pop_county, 
                             t_min, t_max, dt_state_county, prop_gen = 0.6){
  ## Compute the probability of connection between counties using gravity model
  dt_distance[, nb_commut := pop_county1**a*exp(-distance_km*b)]
  dt_distance[distance_km > gamma, nb_commut := 0]
  sum_commut <- dt_distance[, lapply(.SD,sum), by = county2,
                            .SDcols = "nb_commut"]
  setkey(sum_commut, county2)
  dt_distance[, probs := nb_commut / sum_commut[county2,nb_commut]]
  
  # Initialize dt_cases
  dt_cases <- data.table(NULL)
  count <- 1
  i <- 1
  
  while(dim(dt_cases)[1] < nb_cases){
    # Draw the infection date of the ancestors
    date_i <- sample(x = t_min:t_max, size = 1)
    date_i <- as.Date(date_i, origin = "1970-01-01")
    # Draw the county of the ancestors
    county_i <- sample(x =names(pop_county), size = 1, 
                     prob = pop_county/sum(pop_county))
    # Draw the age of the ancestors
    age_group_i <- ceiling(rexp(n = 1, rate = 0.2))
    if(age_group_i>dim(polymod_prop)[2]) age_group_i <- dim(polymod_prop)[2]
    # Draw the genotype of the ancestors
    genotype_i <- sample(c("B3", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "G3", "H1", "H2"),
                         size = 1, 
                         prob = c(0.36, 0.05, 0.10, 0.05, 0.01, 0.01, 0.25, 0.10, 
                                  0.01, 0.05, 0.01))
    # Initialize the cluster number and the generation of the cluster
    clust_i <- i
    gen <- 1
    # Merge all the features and add the ancestor to dt_cases
    case_i <- data.table(as.character(count), dt_state_county[county_i, STNAME],
                         date_i, genotype_i, county_i, age_group_i, TRUE, clust_i,
                         gen)
    colnames(case_i) <- c("ID", "State", "Date", "Genotype", "county", 
                          "age_group", "import", "cluster", "generation")
    dt_cases <- rbind(dt_cases, case_i)
    count <- dim(dt_cases)[1]
    while(count<= dim(dt_cases)[1]){
      # Extract features of index and county
      date_i <- dt_cases[count, Date]
      county_i <- dt_cases[count, county]
      age_group_i <- dt_cases[count, age_group]
      genotype_i <- dt_cases[count, Genotype]
      clust_i <- dt_cases[count, cluster]
      r0_i <- r0_state[dt_state_county[county_i, STNAME]]
      gen_i <- dt_cases[count, generation]
      # Draw number of new cases
      new_cases_i <- rpois(1, r0_i)
      count <- count + 1
      if(new_cases_i > 0){
        # Generate the features of the cases infected by the index
        case_subs <- new_case(sec = new_cases_i, date_index = date_i, 
                              county_index = county_i, 
                              age_group_index = age_group_i, 
                              genotype_index = genotype_i, count = count, 
                              cluster = clust_i, generation = gen_i, 
                              dt_distance = dt_distance,
                              dt_state_county = dt_state_county,
                              polymod_prop = polymod_prop, w = w)
        dt_cases <- rbind(dt_cases, case_subs)
      }
    }
    i <- i + 1 
  }
  # 60% of the genotypes are set to be "not attributed"
  gen_reported <- runif(nrow(dt_cases), 0, 1)
  dt_cases[gen_reported < prop_gen, Genotype := "Not attributed"]
  
  
  dt_distance <- dt_distance[,.(county1, county2, distance_km)]
  setkey(dt_distance, county1)
  distance_matrix <- sapply(unique(dt_distance$county1), 
                            function(X){
    return(dt_distance[county2 == X][, distance_km])
  })
  rownames(distance_matrix) <- colnames(distance_matrix)
  
  vect_pop <- pop_county
  names(vect_pop) <- names(pop_county)
  vect_pop <- vect_pop[rownames(distance_matrix)]
  
  toy_outbreak <- list(cases = dt_cases, distance = distance_matrix,
                        population = vect_pop, age_contact = polymod_prop)
  return(toy_outbreak)
}