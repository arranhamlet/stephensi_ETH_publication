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

run_to_process <- "G0kr0C7AmL"

#After this is done load the models that have run
data_load <- list.files(paste0("output/latin_hypercube/djibouti_fit/", run_to_process, "/"), pattern = "_of_", full.names = T)

all_go <- as.data.frame(rbindlist(sapply(data_load, function(x) read.csv(x, stringsAsFactors = FALSE), simplify = FALSE)))
all_go <- subset(all_go, message %in% c("CONVERGENCE: NORM OF PROJECTED GRADIENT <= PGTOL",
                                        "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"))

#We've taken the 100 best from 200 runs, so top 50%
best_100 <- sort(all_go$likelihood)[1:100]
all_go <- all_go[which(all_go$likelihood %in% best_100), ]

all_go_LHC_subset <- all_go
all_go_LHC_subset$row <- 1:nrow(all_go_LHC_subset)
#Load in the model fits
these_runs <- paste0("output/latin_hypercube/djibouti_fit/", run_to_process, "/LHC_fit/", all_go$row, "_of_", all_go$number_of_rows, "_1.csv")
temp_dependent_EIR <- do.call(rbind, sapply(these_runs, function(x) {
  read.csv(x, stringsAsFactors = FALSE)
}, simplify = FALSE))

all_DJI_outcome_runs <- rbind(temp_dependent_EIR)
row.names(all_DJI_outcome_runs) <- NULL

all_mv_prev_upd <-  all_DJI_outcome_runs %>% 
  group_by(type, time) %>% 
  dplyr::summarise(value_min = quantile(value, 0.025),
                   value_q25 = quantile(value, 0.25),
                   value_q75 = quantile(value, 0.75),
                   value_max = quantile(value, 0.975),
                   value_median = median(value),
                   value = mean(value))

all_mv_prev_upd_df <- as.data.frame(all_mv_prev_upd)

all_mv_prev_upd_df$type <- gsub("mv", "Vector density",
                                gsub("prev", "Prevalence",
                                     all_mv_prev_upd_df$type))

all_mv_prev_upd_df$year <- 2010 + plyr::round_any(all_mv_prev_upd_df$time, 365)/365
all_mv_prev_upd_df$year_dec <- 2010 + all_mv_prev_upd_df$time/365

#Last year values
last_year_values <- all_mv_prev_upd_df[(nrow(all_mv_prev_upd_df) - (365 * 3)): nrow(all_mv_prev_upd_df), ]
last_year_vd <- subset(all_mv_prev_upd_df, type == "Vector density")

last_year_mean <- do.call(rbind, sapply(c("Temperature variation in EIP"), function(x){
  val_name <- names(last_year_vd)[which(grepl("value", names(last_year_vd)))]
  data.frame(type = "Vector density",
             scenario = x,
             quantile = val_name,
             value = sapply(val_name, function(y) mean(last_year_vd[, y])),
             stringsAsFactors = FALSE)
}, simplify = FALSE))

write.csv(last_year_mean, file = "data/malaria_mosquito/DJI_fit_vector_density_ranges.csv", row.names = FALSE)
