## pseudo r2 of beta regressions
## full models
## 9.1.24

# setup -------------------------------------------------------------------
library(tidyverse)

folder_in_m <- "data/processed/betareg_models_v2_relaxed/"
folder_out <- "results/tab_models_v2_relaxed/tab_betareg_r2/"

# load models -------------------------------------------------------------
m_full <- read_rds(str_c(folder_in_m, "models_full.rds"))
mc_full <- read_rds(str_c(folder_in_m, "model_combinations_full.rds"))


# get pseudo r2 -----------------------------------------------------------
model_fit <- mc_full %>% 
  mutate(pseudo.r.squared = map_dbl(m_full, "pseudo.r.squared"),
         pr_nh            = str_c(profile, "_", nat_haz)) %>% 
  select(-c(nat_haz, profile)) %>% 
  pivot_wider(names_from   = "pr_nh",
              names_prefix = "pr2_",
              values_from  = "pseudo.r.squared") %>% 
  select(resp_var_type, simtype, pr2_MP_A, pr2_MP_LED, pr2_IP_A, pr2_IP_LED) %>% 
  arrange(resp_var_type, simtype)

write_delim(model_fit,
            file = str_c(folder_out, "betareg_full_r2.csv"),
            delim = ";")





