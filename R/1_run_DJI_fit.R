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

#By default we use 200
number_of_rows <- 200

all_run_DJI <- expand.grid(number_of_rows = number_of_rows,
                                which_row = 1:number_of_rows,
                                random_string = paste0(stri_rand_strings(1, 10)),
                                seed = 1,
                                prop_cases_stephensi = 1,
                                estimate = "mean", 
                                type = "clinical cases in thousands", #Reference the type column in data/malaria_mosquito/djibouti_malaria_MAP.csv
                                reference = "MAP",
                                population_fit = NA, #If NA this will then take the default values 
                                stringsAsFactors = FALSE)

#Run through each row
#Without high performance cluster support, or some sort of parallelisation this will take a very long time
sapply(1:nrow(all_run_DJI), function(x){
  print(paste0("Run ", x, " of ", number_of_rows))
  do.call(cluster_hypercube_sampling_djibouti, all_run_DJI[x, ])
})


