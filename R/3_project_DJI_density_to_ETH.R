library(odin)
library(readr)
library(odin)
library(stringi)
library(data.table)
library(tidyverse)
library(plyr)
library(raster)
library(lhs)
library(pkgbuild)

#Load in all functions
invisible(sapply(list.files("R/functions/", full.names = T, recursive = T), function(x) source(x)))

#Load in vector density to use
DJI_vd <- read.csv("data/malaria_mosquito/DJI_fit_vector_density_ranges.csv", stringsAsFactors = FALSE)

#These are the pre-existing intervention and EIP values for ETH grouped as per described in the paper
MAP_intervention_values <- read.csv("data/shp/ETH_group_model_parameters.csv")

#Set up runs
#In this example, to save time, we're only going to run for 2 scenarios, the median value and only change net coverage and type
run_scenario <- ETH_projection_run_setup(location_characteristics = MAP_intervention_values[c(12, nrow(MAP_intervention_values)), ],
                                         djibouti_vector_density = DJI_vd,
                                         fit_to = "prevs_2_10",
                                         which_prevalence = "prevs_2_10",
                                         save_location = "output/model_run_density_increase/paper_example_run/",
                                         itn_update = c(0, 0.8), 
                                         irs_update = 0, 
                                         larvicide_update = 0,
                                         after_time = 365 * 25,
                                         net_type = c("Pyrethroid", "PBO-pyrethroid"),
                                         quantile_use = "value_median",
                                         total_LHC_run = 10)

run_scenario$custom_EIR_calibrate <- "data/malaria_mosquito/EIR_LHC_updated_falciparum_set_endemic.csv"
run_scenario$custom_hypercube <- "data/malaria_mosquito/LHC/latin_hypercube_100_updated.csv"

#Run through
sapply(1:nrow(run_scenario), function(x){
  message(paste0("Running ", x, " of ", nrow(run_scenario)))
  do.call(impact_of_introduction_and_intervention, run_scenario[x, ])
})




