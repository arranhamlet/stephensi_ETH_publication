

mle_likelihood_fit <- function(init_density_vec, 
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
                                delayMos,
                                ...){
  
  
  #Now expand to density
  if(any(is.na(init_density_vec))){
    year_value <- as.list(...)
    density_vec <- as.numeric(rep(sapply(paste0("year", 1:nrow(benchmark_data)), function(x) year_value[[which(names(year_value) == x)]]), each = 365))
  } else {
    density_vec <- rep(init_density_vec, each = 365)
  }
  
  # Runs deterministic malaria model
  model_run <- try({
    out <- create_r_model_epidemic(odin_model_path = odin_model_path,
                                   het_brackets = het_brackets,
                                   age = age,
                                   num_int = num_int,
                                   init_ft = init_ft,
                                   init_EIR = init_EIR,
                                   country = NULL,#country,
                                   admin_unit = NULL,#admin_unit,
                                   scalar = scalar,
                                   mu0 = mu0,
                                   Q0 = Q0,
                                   chi = chi,
                                   bites_Bed = bites_Bed,
                                   bites_Indoors = bites_Indoors,
                                   custom_seasonality = custom_seasonality,
                                   time_length = length(density_vec),
                                   density_vec = density_vec,
                                   delayMos = delayMos)
    
    # Extracts the prevelance data
    mod <- out$generator(user = out$state, use_dde = TRUE)
    
    # Runs the model
    mod_run <- mod$run(t = 1:length(density_vec))
    out <- mod$transform_variables(mod_run)
    out_df <- as.data.frame(out)  
    
    # Computes loglikelihood
    year_incidence <- colSums(matrix(out_df$Incidence, ncol = length(density_vec)/365)) * population_size
    LL <- 0
    
    for(i in 1:nrow(benchmark_data)){
      LL <- LL + dnbinom((benchmark_data$number * 1000)[i], 1, 
                         mu = plyr::round_any(year_incidence[i], 1, ceiling), log = TRUE)

    }
    
    
    mod_data <- data.frame(year = 1:length(year_incidence), 
                           type = rep(c("Model", "Data"), each = length(year_incidence)),
                           value = c(year_incidence, benchmark_data$number * 1000))

    message("Likelihood: ", LL)
    
    list(LL, mod_data)
    
  })
  
  if(all(!is.na(init_density_vec))){
    if(class(model_run) == "try-error"){
      list(likelihood = 9^100, model_data_diff = data.frame(year = NA, type = NA, value = NA))
    } else {
      list(likelihood = if(class(model_run) == "try-error") 9^100 else if(MLE_or_MCMC == "MLE") -model_run[[1]] else model_run[[1]],
           model_data_diff = if(class(model_run) == "try-error") NA else model_run[[2]])
    }
  } else {
    if(class(model_run) == "try-error") 9^100 else if(MLE_or_MCMC == "MLE") -model_run[[1]] else model_run[[1]]
  }

}