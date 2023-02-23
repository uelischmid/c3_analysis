## analyze data with boosted regression trees
# full models (across strata)
## 22.2.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)
library(dismo)

folder_in <- "data/processed/analysis/"
folder_out <- "results/brt_models/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in, "analysis_data_transf.rds"))

# LT A
data_LT_A <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "A") %>% 
  as.data.frame()

# LT LED
data_LT_LED <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "LED") %>% 
  as.data.frame()

# ST A
data_ST_A <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "A") %>% 
  as.data.frame()

# ST LED
data_ST_LED <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "LED") %>% 
  as.data.frame()


# build models across strata ----------------------------------------------
model_combinations <- expand_grid(simtype       = c("LT", "ST"),
                                  nat_haz       = c("A", "LED"),
                                  profile       = c("MP", "IP"),
                                  resp_var_type = c("abs", "diff")) %>% 
  mutate(data_in       = str_c("data_", simtype, "_", nat_haz),
         cols_expl_var = case_when(simtype == "LT" ~ "xvar_LT",
                                   TRUE            ~ "xvar_ST"),
         col_resp_var  = case_when(profile == "MP" & resp_var_type == "abs"  ~ 15,
                                   profile == "MP" & resp_var_type == "diff" ~ 16,
                                   profile == "IP" & resp_var_type == "abs"  ~ 17,
                                   profile == "IP" & resp_var_type == "diff" ~ 18,
                                   TRUE ~ 999),
         tree_compl    = 5,
         learning_rate = case_when(simtype == "LT" ~ 0.02,
                                   TRUE            ~ 0.1),
         bag_fraction  = 0.75)

xvar_LT <- c(3, 6, 8, 10, 11, 13)
xvar_ST <- c(3, 6, 8, 9, 10, 11, 13)

models_full <- vector(mode = "list", length = 16)

for (i in seq_along(models_full)) {
  vals <- model_combinations[i, ]
  set.seed(1)
  models_full[[i]] <- gbm.step(data            = get(pull(vals, data_in)),
                               gbm.x           = get(pull(vals, cols_expl_var)),
                               gbm.y           = pull(vals, col_resp_var),
                               family          = "gaussian",
                               tree.complexity = pull(vals, tree_compl),
                               learning.rate   = pull(vals, learning_rate),
                               bag.fraction    = pull(vals, bag_fraction))
}

write_rds(models_full,
          str_c(folder_out, "models_full.rds"))
write_rds(model_combinations,
          str_c(folder_out, "model_combinations_full.rds"))
