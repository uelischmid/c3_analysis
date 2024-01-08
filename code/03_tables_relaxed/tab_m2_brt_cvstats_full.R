## Table of CV-statistics of full BRT models
## 8.1.24


# setup -------------------------------------------------------------------
library(tidyverse)

path_in <- "data/processed/brt_models_v2_relaxed/"
path_out <- "results/tab_models_v2_relaxed/tab_brt_cvstats/"


# load data ---------------------------------------------------------------
brt_full_cv_stats_orig <- read_rds(str_c(path_in, "models_cvstats_full.rds"))


# modify and save table ---------------------------------------------------
brt_full_cv_stats <- brt_full_cv_stats_orig %>% 
  dplyr::select(resp_var_type, simtype, profile, nat_haz, n_trees, learning_rate, deviance.mean, deviance.se, correlation.mean, correlation.se) %>% 
  mutate(profile = factor(profile, levels = c("MP", "IP")),
         nat_haz = factor(nat_haz, levels = c("A", "LED"))) %>% 
  arrange(resp_var_type, simtype, profile, nat_haz)

write_delim(brt_full_cv_stats,
            file = str_c(path_out, "brt_full_cv_stats.csv"),
            delim = ";")
