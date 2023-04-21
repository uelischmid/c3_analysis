
# get raw values ----------------------------------------------------------


library(tidyverse)
# library(betareg)
# library(broom)

folder_in <- "data/processed/nais_analysis_data/"

data_analysis <- read_rds(str_c(folder_in, "analysis_data_transf.rds")) %>% 
  filter(q_reg2 == "normal")

LT_diff_UM <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(stratum == "UM") %>% 
  select(simtype:mgm_intensity_s, sha_i_IP_met_abs, sha_i_IP_met_diff)

LT_diff_UM %>% 
  filter(sha_i_IP_met_diff > 0) %>% 
  mutate(mgm_interval = factor(mgm_interval),
         mgm_intensity = factor(mgm_intensity)) %>% 
  summary()


# get NOM -----------------------------------------------------------------
library(tidyverse)
data_LT_raw <- read_rds("data/raw/nais_analysis_data/LT.rds") %>% 
  select(stratum:nat_haz, starts_with("sha_i"))

LT_NOM_UM_red <- data_LT_raw %>% 
  filter(stratum == "UM") %>% 
  filter(q_reg > 1) %>% 
  filter(mgm == "NOM")
