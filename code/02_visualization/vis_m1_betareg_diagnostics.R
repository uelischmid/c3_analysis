## diagnostic plots of betareg-models
## 22.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v1/"
folder_out <- "results/vis_models_v1/betareg/vis_betareg_diagnostics/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in_d, "analysis_data_transf.rds"))

data_LT_A <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "A")

data_LT_LED <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "LED")

data_ST_A <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "A")

data_ST_LED <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "LED")


# load models -------------------------------------------------------------
models <- read_rds(str_c(folder_in_m, "models.rds"))
model_combinations <- read_rds(str_c(folder_in_m, "model_combinations.rds")) %>% 
  rowid_to_column()



# select models -----------------------------------------------------------
mc_used <- model_combinations %>% 
  filter(link == "logit") %>% 
  filter(formula_nr != 3)


# plot --------------------------------------------------------------------
for (i in 1:nrow(mc_used)) {
  vars <- mc_used[i,]
  
  jpeg(filename = str_c(folder_out,
                        pull(vars, simtype), "_",
                        pull(vars, nat_haz), "_",
                        pull(vars, profile), "_",
                        pull(vars, resp_var_type), "_",
                        "f", pull(vars, formula_nr), ".jpg"),
       width = 1000,
       height = 1000)
  par(mfrow = c(2, 2))
  plot(models[[pull(vars, rowid)]])
  dev.off()
}
