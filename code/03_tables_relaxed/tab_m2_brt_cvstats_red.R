## Table of CV-statistics of reduced BRT models
## 8.11.23


# setup -------------------------------------------------------------------
library(tidyverse)

path_in <- "data/processed/brt_models_v2_relaxed/"
path_out <- "results/tab_models_v2_relaxed/tab_brt_cvstats/"


# load data ---------------------------------------------------------------
LT_red_brt_cvstats <- read_rds(str_c(path_in, "models_cvstats_red_LT.rds")) %>% 
  mutate(init = "LT")
ST_red_brt_cvstats <- read_rds(str_c(path_in, "models_cvstats_red_ST.rds"))


# modify and save table ---------------------------------------------------
brt_red_cv_stats <- bind_rows(ST_red_brt_cvstats, LT_red_brt_cvstats) %>% 
  mutate(init2 = case_when(init == "1" ~ "ST: young",
                           init == "2" ~ "ST: structured",
                           init == "3" ~ "ST: mature",
                           TRUE        ~ "LT")) %>% 
  dplyr::select(resp_var_type, stratum, init2, profile, nat_haz, n_trees, learning_rate, deviance.mean, deviance.se, correlation.mean, correlation.se) %>% 
  mutate(stratum = factor(stratum, levels = c("UM", "HM", "SA")),
         init2 = factor(init2, levels = c("ST: young", "ST: structured", "ST: mature", "LT")),
         profile = factor(profile, levels = c("MP", "IP")),
         nat_haz = factor(nat_haz, levels = c("A", "LED"))) %>% 
  arrange(resp_var_type, stratum, init2, profile, nat_haz)

write_delim(brt_red_cv_stats,
            file = str_c(path_out, "brt_red_cv_stats.csv"),
            delim = ";")
