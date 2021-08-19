
impact_of_introduction_and_intervention <- function(start_prevalence, 
                                                    treatment_cov, 
                                                    delayMos,
                                                    vector_density_increase,
                                                    save_location,
                                                    save_name = NA,
                                                    save_name_base = NA,
                                                    itn_vector_compact,
                                                    irs_vector_compact,
                                                    larvicide_vector_compact,
                                                    t_itn_compact,
                                                    t_irs_compact,
                                                    t_larvicide_compact,
                                                    intervention_decoded,
                                                    surv_bioassay = 0,
                                                    seasonality = F,
                                                    which_prevalence = "prev",
                                                    burn_in_time = 365 * 10,
                                                    introduction_time = 365 * 5,
                                                    after_time = 365 * 10,
                                                    admin_unit = NULL,
                                                    country = NULL,
                                                    
                                                    LHC_on = F,
                                                    LHC_row_run = 1,
                                                    LHC_total_rows = 1,
                                                    seed = 1,
                                                    location = NA, #This is the location as defined in the ethiopia division files - used as a reference
                                                    net_type = "Pyrethroid",
                                                    save_cols = "Sporozoite rate;Incidence;Prevalence;Vector density",
                                                    custom_EIR_calibrate = F,
                                                    custom_hypercube = F){
  
  custom_seasonality <- NA
  
  #Save arguments
  argg <- c(as.list(environment()), list())
  argument_save <- as.data.frame(t(as.data.frame(unlist(argg))))
  
  loop_values <- if(LHC_on == T & grepl(";", LHC_row_run)) as.numeric(strsplit(LHC_row_run, ";")[[1]]) else LHC_row_run
  save_name_all <- paste0(save_name_base, "_loc", location, "_loop", loop_values)
  
  resistance_param <- as.data.frame(fread("data/net_irs/pyrethroid_net_input.csv"))
  resistance_param_net <- subset(resistance_param, resistance_level == round(surv_bioassay, 2) & type == net_type)
  these_EIR_prev_all <- if(custom_EIR_calibrate == F) as.data.frame(fread("output/EIR_LHC_full_run_LHC_giant_v3_1_to_100.csv")) else as.data.frame(fread(custom_EIR_calibrate)) #rbind(EIR_prev_0_to_10, EIR_prev_10_to_20)
  latin_hypercube <- if(custom_hypercube == F) read.csv(paste0("data/malaria_mosquito/LHC/latin_hypercube_", LHC_total_rows, ".csv"), stringsAsFactors = F) else read.csv(custom_hypercube, stringsAsFactors = FALSE)
  
  #Same stuff that happens each run
  #Unpack interventions
  itn_vector <- as.numeric(strsplit(itn_vector_compact, ";")[[1]])
  irs_vector <- as.numeric(strsplit(irs_vector_compact, ";")[[1]])
  larvicide_vector <- as.numeric(strsplit(larvicide_vector_compact, ";")[[1]])
  
  t_vector_itn <- as.numeric(strsplit(t_itn_compact, ";")[[1]])
  t_vector_irs <- as.numeric(strsplit(t_irs_compact, ";")[[1]])
  t_vector_larvicide <- as.numeric(strsplit(t_larvicide_compact, ";")[[1]])
  
  #Work out coverage and population split
  eff_itn_cov <- if(last(itn_vector) == 0) 1e-25 else last(itn_vector)
  eff_irs_cov <- if(last(irs_vector) == 0) 1e-25 else last(irs_vector)
  
  population_split <- rep(NA, 4)
  
  population_split[1] <- (1 - eff_itn_cov) * (1 - eff_irs_cov)  # {No intervention}
  population_split[2] <- eff_itn_cov * (1 - eff_irs_cov) # 	   {ITN only}
  population_split[3] <- (1 - eff_itn_cov) * eff_irs_cov	#      {IRS only}
  population_split[4] <- eff_itn_cov * eff_irs_cov #	   {Both ITN and IRS}
  
  #Which prevalence is closest and lower than the start
  these_EIR_prev <- as.data.frame(these_EIR_prev_all[which(these_EIR_prev_all$delayMos == unique(these_EIR_prev_all$delayMos)[which.min(abs(unique(these_EIR_prev_all$delayMos)) - delayMos)]), ])
  
  this_treatment <- unique(these_EIR_prev$init_ft)[which.min(abs(unique(these_EIR_prev$init_ft) - treatment_cov))]
  
  this_itn <- unique(these_EIR_prev$itn_cov)[which.min(abs(unique(these_EIR_prev$itn_cov) - if(is.na(itn_vector[which(itn_vector != 0)][1])) 0 else itn_vector[2]))]
  this_irs <- unique(these_EIR_prev$irs_cov)[which.min(abs(unique(these_EIR_prev$irs_cov) - if(is.na(irs_vector[which(irs_vector != 0)][1])) 0 else irs_vector[2]))]
  
  #Decide if we want seasonality or not
  if(seasonality != F){
    if(seasonality == "mean"){
      admin_units_seasonal <- readRDS("data/malaria_mosquito/admin_units_seasonal.rds")
      custom_seasonality <- colMeans(seasonal_profile(admin_units_seasonal[which(admin_units_seasonal$country == "Ethiopia"), ]))
      custom_seasonality <- seasonal_profile(as.data.frame(t(colMeans(admin_units_seasonal[which(admin_units_seasonal$country == "Ethiopia"), 4:ncol(admin_units_seasonal)]))))
    } else if(grepl(";", seasonality)){
      seasonality_here <- strsplit(seasonality, ";")[[1]]
      admin_unit <- seasonality_here[1]
      country <- seasonality_here[2]
    }
  }
  
  #Work out start and end prevalence
  use_this_data <- as.data.frame(subset(these_EIR_prev, 
                                        init_ft == this_treatment &
                                          itn_cov == this_itn &
                                          irs_cov == this_irs))
  
  
  #Loop through LHC
  sapply(1:length(loop_values), function(i){
    
    #Wrapped in try because some of the loop values will fail due to parameter combinations that cant be solved - c'est la vie
    try({

      #Save arguments
      save_name <- save_name_all[i]
      argument_save$save_name <- save_name
      argument_save$loop <- loop_values[i]
      
      #Create save location if it doesnt exist alreadear
      if(!dir.exists(paste0(save_location, paste0("model_output/loop_", i, "/")))) dir.create(paste0(save_location, paste0("model_output/loop_", i, "/")), recursive = T)
      if(!dir.exists(paste0(save_location, paste0("function_arguments/loop_", i, "/")))) dir.create(paste0(save_location, paste0("function_arguments/loop_", i, "/")), recursive = T)
      
      write.csv(argument_save, paste0(save_location, "function_arguments/loop_", i, "/function_arguments_", save_name, ".csv"), row.names = FALSE)
      
      #Set up optional LHC
      if(LHC_on == T){
        
        mu0 <- latin_hypercube[loop_values[i], ]$mu0
        Q0 <- latin_hypercube[loop_values[i], ]$Q0
        chi <- latin_hypercube[loop_values[i], ]$chi 
        bites_Bed <- latin_hypercube[loop_values[i], ]$bites_Bed
        bites_Indoors <- latin_hypercube[loop_values[i], ]$bites_Indoors
        
        #Work out how LHC parameters will change the calibration
        use_this_data <- as.data.frame(subset(use_this_data, loop == i))
        
      } else {
        
        mu0 <- 0.12345679
        Q0 <- 0.25
        chi <- 0.5 
        bites_Bed <- 0.4776
        bites_Indoors <- 0.52186
        
      }
      
      #Work out the introduction rate
      this_start_data <- use_this_data[which.min(abs(use_this_data[, which_prevalence] - start_prevalence)), ]
      this_end_data <- use_this_data[which.min(abs(use_this_data$mv0 - (this_start_data$mv0 + vector_density_increase))), ]
      
      values <- sigmoid(seq(-10, 10, length.out = introduction_time))
      
      density_vec <- c(rep(this_start_data$mv0, burn_in_time), 
                       pmin((values * (this_end_data$mv0 - this_start_data$mv0)) + this_start_data$mv0, this_end_data$mv0),
                       rep(this_end_data$mv0, after_time))
      
      #Adjust larviciding because we're targeting different locations for larvae and so
      #it probably wont affect existing anopheles
      density_gain <- max(density_vec) - min(density_vec)
      larvicide_vector_adjusted <- (density_gain * larvicide_vector)/max(density_vec)
      

      #Run the model
      model_run_fit <- create_r_model_epidemic(
        
        #Odin model to use 
        odin_model_path = "odin/odin_model_seasonality_future_intervention.R",
        
        #Number of interventions, heterogeneity and age groups
        num_int = 4,
        het_brackets = 5,
        age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80), 
        
        #Initial prevalence and where it ends converted to density
        init_EIR = this_start_data$EIR,
        density_vec = density_vec,
        time_length = length(density_vec),
        
        #Initial nets, irs and treatment 
        init_ft = this_treatment,
        itn_cov = eff_itn_cov[1],
        irs_cov = eff_irs_cov[1],
        
        #Coverage of interventions
        ITN_IRS_on = t_vector_itn[1],
        IRS_on = if(length(t_vector_irs[which(t_vector_irs != 0)]) == 0) t_vector_irs[which(t_vector_irs != 0)] else 0,
        itn_vector = itn_vector,
        irs_vector = irs_vector,
        larvicide_vector = larvicide_vector_adjusted,
        
        #Timing of interventions
        t_vector_itn = t_vector_itn,
        t_vector_irs = t_vector_irs,
        t_vector_larvicide = t_vector_larvicide,
        
        #Seasonality
        admin_unit = admin_unit,
        country = country,
        custom_seasonality = custom_seasonality,
        
        #Mosquito bionomic guesstimates
        mu0 = mu0, 
        Q0 = Q0, 
        chi = chi, 
        bites_Bed = bites_Bed, 
        bites_Indoors = bites_Indoors,
        delayMos = delayMos,
        
        #Resistance parameters
        d_ITN0 = resistance_param_net$ERG_d_ITN0,
        r_ITN0 = resistance_param_net$ERG_r_ITN0,
        itn_half_life = resistance_param_net$itn_half_life,
        time_resistance = burn_in_time
        
      )
      
      # Edits equilibrium condition with new coverage split. Using the default of 0.25 to all because its the middle ground
      if(any(is.na(population_split))) population_split <- as.numeric(c(0.25, 0.25, 0.25, 0.25))
      mod <- model_run_fit
      
      mod <- edit_equilibrium_varying_nets_irs(mod, population_split)
      
      mod$state$t_vector_irs <- as.numeric(mod$state$t_vector_irs)
      mod$state$t_vector_itn <- as.numeric(mod$state$t_vector_itn)
      mod$state$t_vector_larvicide <- as.numeric(mod$state$t_vector_larvicide)
      
      mod <- mod$generator(user = mod$state, use_dde = TRUE)
      
      mod_run <- mod$run(t = 1:length(density_vec))
      out <- mod$transform_variables(mod_run)
      model_ran <- as.data.frame(out)
      
      #This gets rid of the dead time in the burn in to reduce size
      model_ran_subset <- model_ran
      
      these_columns <- c("Sporozoite rate" = "sporozoite_rate",
                         "Incidence" = "Incidence",
                         "Prevalence" = "prev",
                         "Prevalence_2_10" = "prev_2to10",
                         "Vector density" = "mv")
      
      unpack_take_these <- if(grepl(";", save_cols)) strsplit(save_cols, ";")[[1]] else save_cols
      
      extract_values <- unlist(model_ran_subset[, these_columns[which(names(these_columns) %in% unpack_take_these)]])
      
      name_take <- unique(gsub('[[:digit:]]+', '', names(extract_values)))
      
      total_df <- data.frame(save_name = save_name,
                             location = location,
                             int_decode = intervention_decoded,
                             intervention_combo = paste(c(itn_vector_compact, irs_vector_compact, larvicide_vector_compact, t_itn_compact, t_irs_compact, t_larvicide_compact), collapse = "_"),
                             year = model_ran_subset$t,
                             type = rep(unpack_take_these, 
                                        each = length(model_ran_subset$t)),
                             value = as.numeric(extract_values),
                             stringsAsFactors = FALSE)
      
      fwrite(total_df, paste(save_location, paste0("model_output/loop_", i, "/model_output_", save_name, ".csv"), sep = "/"), row.names = FALSE)
      
    })
  })
}

