## predict PQ values and identify best variable combination
## betareg, reduced models, diff
## 8.11.23


# setup -------------------------------------------------------------------
library(tidyverse)
# library(betareg)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v2_relaxed/"

# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in_d, "analysis_data_transf_relaxed.rds")) %>% 
  filter(q_reg2 == "normal")

# subset data and save as new objects
data_comb <- expand_grid(simtype = "LT",
                         nat_haz = c("A", "LED"),
                         stratum = c("UM", "HM", "SA"))
for (i in 1:nrow(data_comb)) {
  vals <- data_comb[i,]
  object_name <- str_c("data_",
                       pull(vals, simtype), "_",
                       pull(vals, nat_haz), "_",
                       pull(vals, stratum))
  object_content <- data_analysis %>% 
    filter(simtype == pull(vals, simtype)) %>% 
    filter(nat_haz == pull(vals, nat_haz)) %>% 
    filter(stratum == pull(vals, stratum)) %>% 
    mutate(sha_i_MP_met_abs_t = (sha_i_MP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.),
           sha_i_IP_met_abs_t = (sha_i_IP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.))
  
  assign(object_name, object_content)
}

rm(object_name, object_content)


# load models -------------------------------------------------------------
models <- read_rds(str_c(folder_in_m, "models_red_LT.rds"))
model_combinations <- read_rds(str_c(folder_in_m, "model_combinations_red_LT.rds")) %>% 
  rowid_to_column()
model_combinations_diff <- model_combinations %>% 
  filter(resp_var_type == "diff") %>% 
  mutate(init = "LT")


# temp predict ------------------------------------------------------------
fitted_all <- vector(mode = "list", length = nrow(model_combinations_diff))
fitted_best <- fitted_all

for (i in 1:nrow(model_combinations_diff)) {
  mc_curr <- model_combinations_diff[i, ]
  m_curr <- models[[pull(mc_curr, rowid)]]
  
  fitted_all[[i]] <- m_curr$model  %>% 
    as_tibble() %>% 
    mutate(fitted_t  = m_curr$fitted.values,
           fitted_bt = 2 * fitted_t - 1) %>%
    bind_cols(select(mc_curr, -rowid)) %>% 
    select(resp_var_type, stratum, init, profile, nat_haz, q_site2:mgm_intensity, fitted_bt)
  
  fitted_best[[i]] <- fitted_all[[i]] %>% 
    slice_max(fitted_bt, n = 1)
}

fitted_all <- bind_rows(fitted_all)
fitted_best <- bind_rows(fitted_best)
