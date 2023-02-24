## visualize effects of betareg-models
## 20.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(effects)

folder_in_d <- "data/processed/analysis/"
folder_in_m <- "data/processed/betareg_models/"
folder_out <- "results/vis_models_v1/betareg/vis_betareg_effects_singlemodel/"


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
model_combinations <- read_rds(str_c(folder_in_m, "model_combinations.rds"))

# visualize effects -------------------------------------------------------
for (i in seq_along(models)) {
  cat("\r ", i, "/", length(models))
  
  vals <- model_combinations[i, ]
  
  jpeg(filename = str_c(folder_out, "y_auto/br_",
                        pull(vals, simtype), "_",
                        pull(vals, nat_haz), "_",
                        pull(vals, profile), "_",
                        pull(vals, resp_var_type), "_",
                        "f", pull(vals, formula_nr), "_",
                        pull(vals, link), ".jpg"),
       width    = 1000,
       height   = 700)
  plot(allEffects(mod = models[[i]]))
  dev.off()
  
  jpeg(filename = str_c(folder_out, "y_01/br_",
                        pull(vals, simtype), "_",
                        pull(vals, nat_haz), "_",
                        pull(vals, profile), "_",
                        pull(vals, resp_var_type), "_",
                        "f", pull(vals, formula_nr), "_",
                        pull(vals, link), ".jpg"),
       width    = 1000,
       height   = 700)
  plot(allEffects(mod = models[[i]]),
       axes = list(y = list(lim = c(0, 1),
                            type = "response")))
  dev.off()
}