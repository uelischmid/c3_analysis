## diagnostic plots of betareg-models
## reduced models: LT
## relaxed assessments
## 12.5.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v2_relaxed/"
folder_out <- "results/vis_models_v2_relaxed/betareg/vis_betareg_red_diagnostics/"


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


# plot --------------------------------------------------------------------
for (i in 1:nrow(model_combinations)) {
  vals <- model_combinations[i,]
  
  jpeg(filename = str_c(folder_out,
                        pull(vals, simtype), "_",
                        pull(vals, nat_haz), "_",
                        pull(vals, profile), "_",
                        pull(vals, resp_var_type), "_",
                        pull(vals, stratum),".jpg"),
       width = 1000,
       height = 1000)
  par(mfrow = c(2, 2))
  plot(models[[pull(vals, rowid)]], which = 1:4)
  dev.off()
}
