create_r_model_epidemic <- function(odin_model_path,
                                    het_brackets = 5,
                                    age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80), 
                                    country = NULL,
                                    admin_unit = NULL,
                                    time_length,
                                    density_vec,
                                    scalar = 1,
                                    init_ft,
                                    init_EIR,
                                    custom_seasonality = NA,
                                    reduce_betaa = NA,
                                    # time_resistance = burn_in_time,
                                    ...){
  
  ## create model param list using necessary variables
  mpl <- model_param_list_create(...)
  mpl$age_2 <- which(age == 2)
  mpl$age_10 <- which(age == 10)
  
  if(!exists("surv_bioassay")) mpl$surv_bioassay <- 0
  
  # generate initial state variables from equilibrium solution
  state <- epidemic_init_create(age_vector = age, 
                                EIR = init_EIR,
                                ft = init_ft,
                                model_param_list = mpl, 
                                het_brackets = het_brackets,
                                country = NULL,
                                admin_unit = NULL,
                                density_vec = density_vec)
  
  state$density_vec <- as.numeric(density_vec)
  state$reduce_betaa <- if(any(is.na(reduce_betaa))) rep(1, time_length) else as.numeric(reduce_betaa)
  # state$time_resistance <- time_resistance
  # state$d_ITN0_input <- d_ITN0_input
  # state$r_ITN0_input <- r_ITN0_input
  # state$itn_half_life_input <- itn_half_life_input
  
  state$custom_seasonality <- if(all(!is.na(custom_seasonality))){
    
    if(length(rep(1, time_length)) %% 365 == 0){
      rep(custom_seasonality, time_length/365)
    } else {
      c(rep(custom_seasonality, round(time_length/365)), custom_seasonality[1:(length(rep(1, time_length)) %% 365)])
    }
    
  } else if(!is.null(admin_unit)){
    admin_units_seasonal <- readRDS("data/malaria_mosquito/admin_units_seasonal.rds")
    admin_matches <- admin_match(admin_unit = admin_unit, country = country,
                                 admin_units_seasonal = admin_units_seasonal)
    loc_seasonality <- seasonal_profile(admin_units_seasonal[admin_matches, ], scalar = scalar)
    
    if(length(rep(1, time_length)) %% 365 == 0){
      rep(loc_seasonality, time_length/365)
    } else {
      c(rep(loc_seasonality, round(time_length/365)), loc_seasonality[1:(length(rep(1, time_length)) %% 365)])
    }
  } else {
    rep(1, time_length)
  }
  
  mpl$custom_seasonality <- state$custom_seasonality
  
  state$time_length <- time_length
  
  state$age15 <- which(age == 15)
  state$agemax <- which.max(age)
  
  # create odin generator
  gen <- odin::odin(odin_model_path, verbose = FALSE)
  state <- state[which(names(state) %in% coef(gen)[, 1])]
  
  # return mod
  return(list("state" = state, "generator" = gen, "mpl" = mpl))
  
}
