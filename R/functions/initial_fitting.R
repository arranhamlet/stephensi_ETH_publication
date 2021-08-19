
initial_fitting <- function(init_density_vec, 
                                   init_ft,
                                   init_EIR,
                                   
                                   lower_bound = 0.001,
                                   upper_bound = 20,
                                   
                                   country,
                                   admin_unit,
                                   scalar,
                                   mu0,
                                   Q0,
                                   chi,
                                   bites_Bed,
                                   bites_Indoors,
                                   custom_seasonality,
                                   benchmark_data,
                                   population_size,
                                   MLE_or_MCMC,
                                   het_brackets,
                                   age,
                                   num_int,
                                   odin_model_path,
                                   delayMos,
                                   n_tries = 100,
                                   ...){
  
  density_vec_change <- init_density_vec
  likelihood <- rep(NA, n_tries)
  differences <- matrix(rep(NA, n_tries * nrow(benchmark_data)), nrow = n_tries)
  vec_change_history <- matrix(rep(NA, n_tries * nrow(benchmark_data)), nrow = n_tries)
  
  time1 <- Sys.time()
  
  for(i in 1:n_tries){
    
    what_run_go <- mle_likelihood_fit(density_vec_change, 
                                       init_ft,
                                       init_EIR,
                                       country,
                                       admin_unit,
                                       scalar,
                                       mu0,
                                       Q0,
                                       chi,
                                       bites_Bed,
                                       bites_Indoors,
                                       custom_seasonality,
                                       benchmark_data,
                                       population_size,
                                       MLE_or_MCMC,
                                       het_brackets,
                                       age,
                                       num_int,
                                       odin_model_path,
                                       delayMos)
    
    likelihood[i] <- what_run_go[[1]]
    data_diff <- what_run_go[[2]]
    
    if(all(!is.na(data_diff$year))){
      differences[i, ] <- data_diff[which(data_diff$type == "Data"), ]$value - data_diff[which(data_diff$type == "Model"), ]$value
      if(sum(abs(differences[i, ])) < 5000){
        break
      } else if(all(!is.na(data_diff$value))){
        density_vec_change <- density_vec_change + sign(differences[i, ]) * (upper_bound - lower_bound)/50#upper_bound/n_tries
      }
    } else {
      density_vec_change <- vec_change_history[max(which(!is.na(vec_change_history[, 1]))), ] + sign(differences[max(which(!is.na(differences[, 1]))), ]) * sample(seq(0.01, (upper_bound/n_tries) * 2, length.out = 11), 1)
    }
    
    density_vec_change <- pmax(density_vec_change, 0.1)
    vec_change_history[i, ] <- density_vec_change
    
  }
  
  time2 <- Sys.time()
  print(time2 - time1)
  
  list(density_vec_change = density_vec_change,
       differences = differences)
  
}