
run_mle <- function(init_ft,
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
                    init_density_vec = NA,
                    MLE_or_MCMC = "MLE",
                    het_brackets,
                    age,
                    num_int,
                    odin_model_path,
                    ...){
  
  list2env(list(...), environment())
  
  # Start the clock!
  ptm <- proc.time()

  # Sets up the cluster
  lower <- sapply(paste0("lower_year", 1:nrow(benchmark_data)), function(x) get(x))
  upper <- sapply(paste0("upper_year", 1:nrow(benchmark_data)), function(x) get(x))
  
  names(lower) <- gsub("lower_", "", names(lower))
  names(upper) <- gsub("upper_", "", names(upper))
  par <- sapply(paste0("year", 1:nrow(benchmark_data)), function(x) get(x))
  
  message("Starting initial fitting")
  
  #Run initial shoddy fitting to get closer to the real number
  par_updated <- initial_fitting(init_density_vec = par, 
                                        init_ft,
                                        init_EIR,
                                        lower_bound = 0.001,
                                        upper_bound = max(upper),
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
                                        n_tries = 10)
  
  message("Starting actual fitting")
  
  #Run
  names(par_updated[[1]]) <- paste0("year", 1:length(par_updated[[1]]))
  
  theta_parallel <- optim(par = par_updated[[1]], 
                          fn = mle_likelihood_fit, 
                          gr = NULL,
                          method  = "BFGS",#"L-BFGS-B",
                          control = list(factr = 1e12),
                          lower = lower,
                          upper = upper, 
                          init_density_vec = NA,
                          init_ft = init_ft,
                          init_EIR = init_EIR,
                          country = country,
                          admin_unit = admin_unit,
                          scalar = scalar,
                          mu0 = mu0,
                          Q0 = Q0,
                          chi = chi,
                          bites_Bed = bites_Bed,
                          bites_Indoors = bites_Indoors,
                          custom_seasonality = custom_seasonality,
                          benchmark_data = benchmark_data,
                          population_size = population_size,
                          MLE_or_MCMC = MLE_or_MCMC,
                          het_brackets = het_brackets,
                          age = age,
                          num_int = num_int,
                          odin_model_path = odin_model_path,
                          delayMos = delayMos)
  
  # Stop the clock
  print(proc.time() - ptm)
  
  list(theta_parallel, differences = par_updated[[2]])
  
}