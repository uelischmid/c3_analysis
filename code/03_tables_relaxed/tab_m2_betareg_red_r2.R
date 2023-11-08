## pseudo r2 of beta regressions
## reduced models
## 8.11.23

# setup -------------------------------------------------------------------
library(tidyverse)

folder_in_m <- "data/processed/betareg_models_v2_relaxed/"
folder_out <- "results/tab_models_v2_relaxed/tab_betareg_r2/"

# load models -------------------------------------------------------------
m_ST <- read_rds(str_c(folder_in_m, "models_red_ST.rds"))
mc_ST <- read_rds(str_c(folder_in_m, "model_combinations_red_ST.rds")) %>% 
  rowid_to_column()

m_LT <- read_rds(str_c(folder_in_m, "models_red_LT.rds"))
mc_LT <- read_rds(str_c(folder_in_m, "model_combinations_red_LT.rds")) %>% 
  rowid_to_column() %>% 
  mutate(init = "LT")


# get pseudo r2 -----------------------------------------------------------
model_fit <- bind_rows(mc_ST %>%
                         mutate(pseudo.r.squared = map_dbl(m_ST, "pseudo.r.squared")),
                       mc_LT %>% 
                         mutate(pseudo.r.squared = map_dbl(m_LT, "pseudo.r.squared"))) %>% 
  mutate(init2 = case_when(init == "1" ~ "ST: young",
                           init == "2" ~ "ST: structured",
                           init == "3" ~ "ST: mature",
                           TRUE        ~ "LT"),
         init2 = factor(init2, levels = c("ST: young", "ST: structured", "ST: mature", "LT")),
         pr_nh = str_c(profile, "_", nat_haz),
         stratum = factor(stratum, levels = c("UM", "HM", "SA"))) %>% 
  select(-c(rowid, simtype, nat_haz, profile, init)) %>% 
  pivot_wider(names_from = "pr_nh",
              names_prefix = "pr2_",
              values_from = "pseudo.r.squared") %>% 
  arrange(resp_var_type, stratum, init2) %>% 
  select(resp_var_type, stratum, init2, pr2_MP_A, pr2_MP_LED, pr2_IP_A, pr2_IP_LED)


write_delim(model_fit,
            file = str_c(folder_out, "betareg_red_r2.csv"),
            delim = ";")





