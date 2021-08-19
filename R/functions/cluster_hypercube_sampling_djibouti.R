
cluster_hypercube_sampling_djibouti <- function(number_of_rows, which_row, random_string, seed, prop_cases_stephensi = 1,
                                                estimate = "mean", type = "clinical cases in thousands", reference = "MAP",
                                                population_fit = NA, custom_seasonality = NA){
  
  #Load in data and set up parameters
  dji_malaria_data <- read.csv("data/malaria_mosquito/djibouti_malaria_MAP.csv", stringsAsFactors = FALSE)
  dji_population_data <- read.csv("data/malaria_mosquito/DJI_population_at_risk.csv", stringsAsFactors = FALSE)
  
  #Load in EIR data
  EIR_prev <- as.data.frame(fread("data/malaria_mosquito/relationship_of_EIR_intervention_model_output.csv", stringsAsFactors = FALSE))
  
  #Create LHC
  all_values <- list(
    mu0 = rnorm(number_of_rows, mean = 0.12345679, sd = 0.12345679 * 0.25), #Daily mortality
    Q0 = seq(.1, .4, length.out = number_of_rows), #Anthropophagy
    chi = rnorm(number_of_rows, mean = 0.5, sd = 0.5 * 0.25), #Endophily
    bites_Indoors = rnorm(number_of_rows, mean = 0.52186, sd = 0.52186 * 0.25), #Proportion of bites indoors
    bites_Bed = rnorm(number_of_rows, mean = 0.4776, sd = 0.4776 * 0.25), #Proportion of bites indoors that are in bed
    delayMos = 10 #Extrinsic incubation period
  )
  
  latin_hypercube <- as.data.frame(randomLHS(length(all_values[[1]]), length(all_values)))
  
  for(i in 1:length(all_values)){
    latin_hypercube[, i] <- qunif(latin_hypercube[, i], min(all_values[[i]]), max(all_values[[i]]))
  }
  
  colnames(latin_hypercube) <- names(all_values)
  
  #Create directory for the save files
  if(!dir.exists(paste0("output/latin_hypercube/djibouti_fit/", random_string, "/"))) dir.create(paste0("output/latin_hypercube/djibouti_fit/", random_string, "/"), recursive = T)
  
  #Save LHC for reference
  write.csv(latin_hypercube, 
            file = paste0("output/latin_hypercube/djibouti_fit/", random_string, "/index_latin_hypercube_", random_string, ".csv"), row.names = FALSE)
  
  #Find the closest EIR to the prevalence for that level of anthropophagy
  these_EIR_prev <- EIR_prev[which(EIR_prev$delayMos == latin_hypercube[which_row, ]$delayMos), ]
  these_EIR_prev <- these_EIR_prev[which.min(these_EIR_prev$Q0 - latin_hypercube[which_row, ]$Q0), ]
  
  #Work out the number of cases that we are attributing to stephensi - default is 1
  dji_malaria_data$number <- dji_malaria_data$number * prop_cases_stephensi 
  
  #Run the model fitting
  run_mle_go <- run_mle(
    
    #Model dimensions
    num_int = 4,
    het_brackets = 5,
    age = c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3.5, 5, 7.5, 10, 15, 20, 30, 40, 50, 60, 70, 80),
    
    #Model parameters to work out transmission
    init_ft = 0,
    itn_cov = 0.1,
    irs_cov = 0,
    
    #Input parameters
    delayMos = latin_hypercube[which_row, ]$delayMos,
    
    #Set up country start EIR
    init_EIR = these_EIR_prev[which.min(abs(these_EIR_prev$prev - 0.005)), ]$EIR,
    country = "Djibouti",
    admin_unit = "Djibouti",
    scalar = 1,
    
    #LHC bionomics
    mu0 = latin_hypercube[which_row, ]$mu0, 
    Q0 = latin_hypercube[which_row, ]$Q0, 
    chi = latin_hypercube[which_row, ]$chi, 
    bites_Bed = latin_hypercube[which_row, ]$bites_Bed, 
    bites_Indoors = latin_hypercube[which_row, ]$bites_Indoors, 
    custom_seasonality = custom_seasonality,
    
    #Benchmark data to compare model outputs to
    benchmark_data = dji_malaria_data[which(dji_malaria_data$estimate == estimate &
                                              dji_malaria_data$type == type &
                                              dji_malaria_data$reference == reference), ],
    population_size = if(is.na(population_fit)) subset(dji_population_data, year %in% unique(dji_malaria_data$year))$population_at_risk else 1000000,
    init_density_vec = NA,
    MLE_or_MCMC = "MLE",
    
    #Model output
    odin_model_path = "odin/odin_model_seasonality.R",
    
    #Now guess the vector density and bounds to change to fit to the incidence
    year1 = 1,
    year2 = 1.25,
    year3 = 1.5,
    year4 = 8,
    year5 = 7,
    year6 = 4,
    year7 = 3.5,
    year8 = 3.25,
    year9 = 3.75,
    year10 = 4.25,
    
    lower_year1 = 0.001,
    lower_year2 = 0.001,
    lower_year3 = 0.001,
    lower_year4 = 0.001,
    lower_year5 = 0.001,
    lower_year6 = 0.001,
    lower_year7 = 0.001,
    lower_year8 = 0.001,
    lower_year9 = 0.001,
    lower_year10 = 0.001,
    
    upper_year1 = 90,
    upper_year2 = 90,
    upper_year3 = 90,
    upper_year4 = 90,
    upper_year5 = 90,
    upper_year6 = 90,
    upper_year7 = 90,
    upper_year8 = 90,
    upper_year9 = 90,
    upper_year10 = 90
    
  )
  
  message("Saving output")
  
  #Save output
  differences <- t(data.frame(run_mle_go[[2]][nrow(run_mle_go[[2]]), ]))
  colnames(differences) <- paste0("data_model_difference_", 1:ncol(differences))
  row.names(differences) <- NULL
  
  save_this <- cbind(data.frame(number_of_rows = number_of_rows, 
                                row = which_row,
                                prop_cases_stephensi =prop_cases_stephensi,
                                message = run_mle_go[[1]]$message,
                                likelihood = run_mle_go[[1]]$value), 
                     as.data.frame(latin_hypercube[which_row, ]),
                     differences,
                     t(as.data.frame(run_mle_go[[1]]$par)))
  
  #Write output
  write.csv(save_this, 
            file = paste0("output/latin_hypercube/djibouti_fit/", random_string, "/", which_row, "_of_", number_of_rows, "_", random_string, "_", prop_cases_stephensi,".csv"), 
            row.names = FALSE)
  
  
}