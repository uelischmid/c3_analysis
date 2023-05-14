## analyze data with boosted regression trees
## models with reduced data: ST
## relaxed assessments
## 12.5.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)
library(dismo)

folder_in <- "data/processed/nais_analysis_data/"
folder_out <- "data/processed/brt_models_v2_relaxed/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in, "analysis_data_transf_relaxed.rds")) %>% 
  filter(q_reg2 == "normal")

data_comb <- expand_grid(simtype = "ST",
                         nat_haz = c("A", "LED"),
                         stratum = c("UM", "HM", "SA"),
                         init    = c("1", "2", "3"))
for (i in 1:nrow(data_comb)) {
  vals <- data_comb[i,]
  object_name <- str_c("data_",
                       pull(vals, simtype), "_",
                       pull(vals, nat_haz), "_",
                       pull(vals, stratum), "_",
                       pull(vals, init))
  object_content <- data_analysis %>% 
    filter(simtype == pull(vals, simtype)) %>% 
    filter(nat_haz == pull(vals, nat_haz)) %>% 
    filter(stratum == pull(vals, stratum)) %>% 
    filter(init == pull(vals, init)) %>% 
    as.data.frame()
  
  assign(object_name, object_content)
}

rm(object_name, object_content)


# build models across strata ----------------------------------------------
model_combinations <- expand_grid(simtype       = c("ST"),
                                  nat_haz       = c("A", "LED"),
                                  profile       = c("MP", "IP"),
                                  resp_var_type = c("abs", "diff"),
                                  stratum       = c("UM", "HM", "SA"),
                                  init          = c("1", "2", "3")) %>% 
  mutate(data_in       = str_c("data_", simtype, "_", nat_haz, "_", stratum, "_", init),
         cols_expl_var = "xvar",
         col_resp_var  = case_when(profile == "MP" & resp_var_type == "abs"  ~ 15,
                                   profile == "MP" & resp_var_type == "diff" ~ 16,
                                   profile == "IP" & resp_var_type == "abs"  ~ 17,
                                   profile == "IP" & resp_var_type == "diff" ~ 18,
                                   TRUE ~ 999),
         tree_compl    = 3,
         learning_rate = c(0.01, 0.01, 0.01, 0.1, 0.007, 0.02, 0.01, 0.01, 0.01, # ST A MP abs UM (1, 2, 3) - SA (1, 2, 3) ## CHECK AGAIN
                           0.01, 0.01, 0.01, 0.1, 0.007, 0.01, 0.007, 0.01, 0.01, # ST A MP diff
                           0.008, 0.009, 0.1, 0.01, 0.03, 0.01, 0.1, 0.01, 0.01, # ST A IP abs
                           0.007, 0.007, 0.1, 0.01, 0.03, 0.01, 0.01, 0.01, 0.01, # ST A IP diff
                           0.005, 0.005, 0.01, 0.01, 0.008, 0.1, 0.005, 0.007, 0.007, # ST LED MP abs
                           0.005, 0.01, 0.01, 0.01, 0.01, 0.1, 0.005, 0.01, 0.01, # ST LED MP diff
                           0.007, 0.01, 0.007, 0.01, 0.01, 0.01, 0.02, 0.01, 0.01, # ST LED IP abs
                           0.01, 0.01, 0.007, 0.008, 0.01, 0.01, 0.02, 0.01, 0.009), # ST LED IP diff
         bag_fraction  = 0.75)

xvar <- c(6, 10, 11, 13)

models <- vector(mode = "list", length = nrow(model_combinations))
models_cvstats <- models

for (i in seq_along(models)) {
  # i <- i + 1
  vals <- model_combinations[i, ]
  set.seed(1)
  mod <- gbm.step(data            = get(pull(vals, data_in)),
                  gbm.x           = get(pull(vals, cols_expl_var)),
                  gbm.y           = pull(vals, col_resp_var),
                  family          = "gaussian",
                  tree.complexity = pull(vals, tree_compl),
                  learning.rate   = pull(vals, learning_rate),
                  # learning.rate   = 0.01,
                  bag.fraction    = pull(vals, bag_fraction))
  
  attr(mod$Terms, ".Environment") <- NULL
  
  models[[i]] <- mod
  
  models_cvstats[[i]] <- vals %>% 
    # dplyr::select(simtype:init) %>% 
    bind_cols(bind_cols(mod$cv.statistics)) %>% 
    mutate(n_trees = length(mod$trees))
}

map(models, "n.trees")

models_cvstats <- bind_rows(models_cvstats)

summary(models_cvstats)

ggplot(models_cvstats, aes(simtype, deviance.mean)) +
  geom_boxplot() +
  facet_wrap(~resp_var_type)

ggplot(models_cvstats, aes(simtype, correlation.mean)) +
  geom_boxplot() +
  facet_wrap(~resp_var_type)

write_rds(models,
          str_c(folder_out, "models_red_ST.rds"))
write_rds(models_cvstats,
          str_c(folder_out, "models_cvstats_red_ST.rds"))
write_rds(model_combinations,
          str_c(folder_out, "model_combinations_red_ST.rds"))


