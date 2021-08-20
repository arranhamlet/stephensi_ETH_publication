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
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

#Read in key
run_key <- read.csv(max(list.files("output/model_run_density_increase/paper_example_run/function_arguments/", ".csv", full.names = T)), stringsAsFactors = FALSE)

#Read in all model runs
all_model_output <- list.files("output/model_run_density_increase/paper_example_run/model_output/", full.names = T, recursive = T)
all_model_run <- as.data.frame(rbindlist(sapply(all_model_output, function(x){
  df <- read.csv(x)
  df$loop <- last(strsplit(df$save_name[1], "loop")[[1]])
  df
}, simplify = FALSE)))

all_model_run$year_relative <- all_model_run$year - run_key$burn_in_time[1]/365

#If we aggregate over locations, we combine LHC values to get an overall value with confidence intervals for each location
#If we aggregate over LHC values, we combine locations to get an overall value with confidence intervals for each LHC draw
results_by_location <- all_model_run %>%
  group_by(location, int_decode, year_relative, type) %>% 
  dplyr::summarise(x = quantile(value, c(0.025, 0.5, 0.975)), q = c("min", "median", "max"))

results_by_location_wide <- as.data.frame(spread(results_by_location, q, x))

results_by_loop <- all_model_run %>%
  group_by(loop, int_decode, type) %>% 
  dplyr::summarise(x = quantile(value, c(0.025, 0.5, 0.975)), q = c("min", "median", "max"))

results_by_loop_wide <- as.data.frame(spread(results_by_loop, q, x))

#Plot the results
incidence_per_hundred_thousand <- ggplot(data = subset(results_by_location_wide, type == "Incidence"), 
       aes(x = year_relative, y = median * 100000, ymin = min * 100000, ymax = max * 100000, fill = int_decode,
           color = int_decode, group = int_decode)) +
  geom_line(size = 2) +
  geom_ribbon(alpha = 0.25) +
  facet_wrap(~location) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Year", y = "Value", fill = "Intervention", color = "Intervention") +
  geom_vline(xintercept = c(0, run_key$introduction_time[1]/365), linetype = "dashed")

prevalence <- ggplot(data = subset(results_by_location_wide, type == "Incidence"), 
                                         aes(x = year_relative, y = median * 100000, ymin = min * 100000, ymax = max * 100000, fill = int_decode,
                                             color = int_decode, group = int_decode)) +
  geom_line(size = 2) +
  geom_ribbon(alpha = 0.25) +
  facet_wrap(~location) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  labs(x = "Year", y = "Value", fill = "Intervention", color = "Intervention") +
  geom_vline(xintercept = c(0, run_key$introduction_time[1]/365), linetype = "dashed")

ggarrange(prevalence, incidence_per_hundred_thousand, ncol = 1, common.legend = T, legend = "bottom")


#Plot the overall results
ggplot(data = subset(results_by_loop_wide, type == "Incidence")) +
  geom_boxplot(aes(x = int_decode, y = median * 100000, fill = type), width = 0.1) +
  geom_violin(aes(x = int_decode, y = median * 100000, fill = type), alpha = 0.5) +
  geom_jitter(aes(x = int_decode, y = median * 100000), height = 0, width = 0.1) +
  facet_wrap(~type, scales = "free_y") +
  theme_minimal() +
  labs(x = "", y = "Median value") +
  theme(legend.position = "none")











