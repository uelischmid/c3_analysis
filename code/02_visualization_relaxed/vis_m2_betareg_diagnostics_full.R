## diagnostic plots of betareg-models
## full models
## relaxed assessments
## 12.5.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v2_relaxed/"
folder_out <- "results/vis_models_v2_relaxed/betareg/vis_betareg_full_diagnostics/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in_d, "analysis_data_transf_relaxed.rds"))

data_LT_A <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "A") %>% 
  mutate(sha_i_MP_met_abs_t = (sha_i_MP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.),
         sha_i_IP_met_abs_t = (sha_i_IP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.))

data_LT_LED <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "LED") %>% 
  mutate(sha_i_MP_met_abs_t = (sha_i_MP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.),
         sha_i_IP_met_abs_t = (sha_i_IP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.))

data_ST_A <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "A") %>% 
  mutate(sha_i_MP_met_abs_t = (sha_i_MP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.),
         sha_i_IP_met_abs_t = (sha_i_IP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.))

data_ST_LED <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "LED") %>% 
  mutate(sha_i_MP_met_abs_t = (sha_i_MP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.),
         sha_i_IP_met_abs_t = (sha_i_IP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.))


# load models -------------------------------------------------------------
models <- read_rds(str_c(folder_in_m, "models_full.rds"))
model_combinations <- read_rds(str_c(folder_in_m, "model_combinations_full.rds")) %>% 
  rowid_to_column()


# plot --------------------------------------------------------------------
for (i in 1:nrow(model_combinations)) {
  vars <- model_combinations[i,]
  
  jpeg(filename = str_c(folder_out,
                        pull(vars, simtype), "_",
                        pull(vars, nat_haz), "_",
                        pull(vars, profile), "_",
                        pull(vars, resp_var_type), ".jpg"),
       width = 1000,
       height = 1000)
  par(mfrow = c(2, 2))
  plot(models[[pull(vars, rowid)]], which = 1:4)
  dev.off()
}
