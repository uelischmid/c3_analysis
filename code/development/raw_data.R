
# get raw values ----------------------------------------------------------
library(tidyverse)
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

LT_SA_IP_abs <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(stratum == "SA") %>% 
  select(nat_haz, q_site2, mgm_type, mgm_interval, mgm_intensity, sha_i_IP_met_abs) %>% 
  pivot_wider(names_from = nat_haz, values_from = sha_i_IP_met_abs) %>% 
  mutate(diff_A_LED = A - LED)

LT_SA_IP_diff <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(stratum == "SA") %>% 
  select(nat_haz, q_site2, mgm_type, mgm_interval, mgm_intensity, sha_i_IP_met_diff) %>% 
  pivot_wider(names_from = nat_haz, values_from = sha_i_IP_met_diff) %>% 
  mutate(diff_A_LED = A - LED)


# get NOM -----------------------------------------------------------------
library(tidyverse)
simtype <- "LT"
stra <- "SA"

data_raw <- read_rds(str_c("data/raw/nais_analysis_data/", simtype, ".rds")) %>% 
  select(stratum:nat_haz, starts_with("sha_i"))

NOM_red <- data_raw %>% 
  filter(stratum == stra) %>% 
  filter(q_reg > 1) %>% 
  filter(mgm == "NOM")
