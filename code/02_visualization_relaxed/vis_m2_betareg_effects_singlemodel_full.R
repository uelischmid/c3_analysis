## visualize effects of betareg full models v2
## relaxed assessments
## 12.5.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(ggeffects)
library(cowplot)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v2_relaxed/"
folder_out <- "results/vis_models_v2_relaxed/betareg/vis_betareg_full_effects_singlemodel/"


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


# plot effects ------------------------------------------------------------
for (i in 1:nrow(model_combinations)) {
  vals <- model_combinations[i,]
  mod <- models[[i]]
  
  effs <- ggeffect(model = mod)
  eff_plots <- plot(effs)
  comb_plot <- plot_grid(plotlist = eff_plots)
  ggsave(str_c(folder_out,
               pull(vals, simtype), "_",
               pull(vals, nat_haz), "_",
               pull(vals, profile), "_",
               pull(vals, resp_var_type), ".jpg"),
         scale = 1.5)
}
