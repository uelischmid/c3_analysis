## analyze data with boosted regression trees
# full models (across strata)
## 22.2.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)
library(dismo)

folder_in <- "data/processed/nais_analysis_data/"
folder_out <- "data/processed/brt_models_v2/"


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
         learning_rate = c(0.05, 0.06, 0.03, 0.02,
                           0.04, 0.04, 0.03, 0.02,
                           0.1, 0.1, 0.1, 0.1,
                           0.1, 0.1, 0.1, 0.1),
         
         bag_fraction  = 0.75)

xvar_LT <- c(3, 6, 8, 10, 11, 13)
xvar_ST <- c(3, 6, 8, 9, 10, 11, 13)

models <- vector(mode = "list", length = nrow(model_combinations))
models_cvstats <- models

for (i in seq_along(models)) {
  vals <- model_combinations[i, ]
  set.seed(1)
  mod <- gbm.step(data            = get(pull(vals, data_in)),
                  gbm.x           = get(pull(vals, cols_expl_var)),
                  gbm.y           = pull(vals, col_resp_var),
                  family          = "gaussian",
                  tree.complexity = pull(vals, tree_compl),
                  learning.rate   = pull(vals, learning_rate),
                  bag.fraction    = pull(vals, bag_fraction))
  
  attr(mod$Terms, ".Environment") <- NULL
  
  models[[i]] <- mod
  
  models_cvstats[[i]] <- vals %>% 
    dplyr::select(simtype:resp_var_type) %>% 
    bind_cols(bind_cols(mod$cv.statistics))
}

models_cvstats <- bind_rows(models_cvstats)

summary(models_cvstats)

ggplot(models_cvstats, aes(simtype, deviance.mean)) +
  geom_boxplot() +
  facet_wrap(~resp_var_type)

ggplot(models_cvstats, aes(simtype, correlation.mean)) +
  geom_boxplot() +
  facet_wrap(~resp_var_type)

write_rds(models,
          str_c(folder_out, "models_full.rds"))
write_rds(models_cvstats,
          str_c(folder_out, "models_cvstats_full.rds"))
write_rds(model_combinations,
          str_c(folder_out, "model_combinations_full.rds"))


