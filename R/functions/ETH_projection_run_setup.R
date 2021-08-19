
ETH_projection_run_setup <- function(location_characteristics,
                              djibouti_vector_density,
                              fit_to,
                              which_prevalence,
                              save_location,
                              resistance = 0.4314444,
                              itn_update = c(0), 
                              irs_update = c(0), 
                              larvicide_update = c(0), 
                              burn_in_time = 365 * 12, 
                              introduction_time = 365 * 3, 
                              after_time = 365 * 15,
                              LHC_on = T,
                              total_LHC_run = 100,
                              LHC_row_run = NA,
                              t_itn_compact = NA,
                              t_irs_compact = NA,
                              t_larvicide_compact = NA,
                              preload_calibrate = NA,
                              net_type = "Pyrethroid",
                              save_cols = "Incidence;Prevalence",
                              quantile_use = c("value_median", "value_min", "value_max"),
                              seasonality = F){
  
  LHC_these <- if(is.na(LHC_row_run)) paste(1:total_LHC_run, collapse = ";") else LHC_row_run
  
  all_location_run <- do.call(rbind, sapply(1:nrow(location_characteristics), function(x){
    
    itn_go <- unique(paste0("0;", location_characteristics$itn[x], ";", pmax(location_characteristics$itn[x], itn_update))) 
    irs_go <- unique(paste0("0;", location_characteristics$irs[x], ";", pmax(location_characteristics$irs[x], irs_update)))
    larvicide_go <- unique(paste0("0;0;", larvicide_update))
    
    #Set up all runs
    values_go <- expand.grid(start_prevalence = location_characteristics[, colnames(location_characteristics) == fit_to][x],
                             treatment_cov = location_characteristics$treatment[x],
                             delayMos = location_characteristics$EIP[x],
                             vector_density_increase = subset(djibouti_vector_density, scenario == "Temperature variation in EIP" & quantile %in% quantile_use)$value,
                             save_location = save_location,
                             which_prevalence = which_prevalence,
                             
                             itn_vector_compact = itn_go,
                             irs_vector_compact = irs_go,
                             larvicide_vector_compact = larvicide_go,
                             
                             t_itn_compact = if(any(is.na(t_itn_compact))) paste0("-100;2;", burn_in_time + introduction_time) else t_itn_compact,
                             t_irs_compact = if(any(is.na(t_irs_compact))) paste0("-100;2;", burn_in_time + introduction_time) else t_irs_compact,
                             t_larvicide_compact = if(any(is.na(t_larvicide_compact))) paste0("-100;2;", burn_in_time + introduction_time) else t_larvicide_compact,
                             
                             surv_bioassay = resistance,
                             
                             seasonality = seasonality,
                             
                             burn_in_time = burn_in_time,
                             introduction_time = introduction_time,
                             after_time = after_time,
                             
                             LHC_on = LHC_on,
                             LHC_row_run = LHC_these,
                             LHC_total_rows = total_LHC_run,
                             seed = 1,
                             location = location_characteristics$unique_id[x],
                             net_type = net_type,
                             save_cols = save_cols,
                             stringsAsFactors = FALSE)
    
    #Decode intervention
    itn_these <- matrix(unlist(strsplit(values_go$itn_vector_compact, ";")), ncol = 3, byrow = T)
    change_itn <- itn_these[, 2] == itn_these[, 3]
    itn_decoded <- change_itn
    itn_decoded[which(change_itn == TRUE)] <- 0
    itn_decoded[which(change_itn == FALSE)] <- as.numeric(itn_these[which(change_itn == FALSE), 3])
    
    irs_these <- matrix(unlist(strsplit(values_go$irs_vector_compact, ";")), ncol = 3, byrow = T)
    change_irs <- irs_these[, 2] == irs_these[, 3]
    irs_decoded <- change_irs
    irs_decoded[which(change_irs == TRUE)] <- 0
    irs_decoded[which(change_irs == FALSE)] <- as.numeric(irs_these[which(change_irs == FALSE), 3])
    
    values_go$intervention_decoded <- paste0("ITN ", itn_decoded * 100, "%/",
                                             "IRS ", irs_decoded * 100, "%/",
                                             "Larvicide ", as.numeric(matrix(unlist(strsplit(values_go$larvicide_vector_compact, ";")), ncol = 3, byrow = T)[, 3])  * 100, "%")
    
    values_go
    
  }, simplify = FALSE))
  
  message(nrow(all_location_run))
  all_location_run$save_name_base <- stri_rand_strings(n = nrow(all_location_run), length = 10)
  if(!dir.exists(paste0(all_location_run$save_location[1], "function_arguments/"))) dir.create(paste0(all_location_run$save_location[1], "function_arguments/"), recursive = T)
  write.csv(all_location_run, paste0(all_location_run$save_location[1], paste0("function_arguments/key_all_location_run_", gsub("-| |:", "", Sys.time())), ".csv"), row.names = FALSE)
  all_location_run
  
}
