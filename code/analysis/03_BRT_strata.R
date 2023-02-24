## analyze data with boosted regression trees
# full models (across strata)
## 22.2.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)
library(dismo)

folder_in <- "data/processed/analysis/"
folder_out <- "data/processed/brt_models/"


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


# build stratum-specific models -------------------------------------------
model_combinations <- expand_grid(simtype       = c("LT", "ST"),
                                  nat_haz       = c("A", "LED"),
                                  profile       = c("MP", "IP"),
                                  resp_var_type = c("abs", "diff"),
                                  stratum       = c("UM", "HM", "SA")) %>% 
  mutate(data_in_part  = str_c("data_", simtype, "_", nat_haz),
         cols_expl_var = case_when(simtype == "LT" ~ "xvar_LT",
                                   TRUE            ~ "xvar_ST"),
         col_resp_var  = case_when(profile == "MP" & resp_var_type == "abs"  ~ 15,
                                   profile == "MP" & resp_var_type == "diff" ~ 16,
                                   profile == "IP" & resp_var_type == "abs"  ~ 17,
                                   profile == "IP" & resp_var_type == "diff" ~ 18,
                                   TRUE ~ 999),
         tree_compl    = 4,
         learning_rate = c(0.010, 0.010, 0.010, # LT A MP abs: UM, HM, SA
                           0.010, 0.010, 0.010, # LT A MP diff
                           0.003, 0.003, 0.003, # LT A IP abs
                           0.003, 0.003, 0.003, # LT A IP diff
                           0.010, 0.010, 0.010, # LT LED
                           0.010, 0.010, 0.010,
                           0.003, 0.003, 0.003,
                           0.003, 0.003, 0.003,
                           0.100, 0.070, 0.090, # ST A
                           0.100, 0.100, 0.100,
                           0.100, 0.100, 0.060,
                           0.100, 0.100, 0.050,
                           0.080, 0.080, 0.060, # ST LED
                           0.090, 0.100, 0.080,
                           0.080, 0.060, 0.080,
                           0.080, 0.080, 0.080),
         bag_fraction  = 0.75)

xvar_LT <- c(6, 8, 10, 11, 13)
xvar_ST <- c(6, 8, 9, 10, 11, 13)


models_stratum <- vector(mode = "list", length = nrow(model_combinations))

for (i in seq_along(models_stratum)) {
  vals <- model_combinations[i,]
  data_in <- get(pull(vals, data_in_part)) %>% 
    filter(stratum == pull(vals, stratum))
  set.seed(1)
  models_stratum[[i]] <- gbm.step(data            = data_in,
                 gbm.x           = get(pull(vals, cols_expl_var)),
                 gbm.y           = pull(vals, col_resp_var),
                 family          = "gaussian",
                 tree.complexity = pull(vals, tree_compl),
                 learning.rate   = pull(vals, learning_rate),
                 bag.fraction    = pull(vals, bag_fraction))
}

for (i in seq_along(models_stratum)) {
  attr(models_stratum[[i]]$Terms, ".Environment") <- NULL
}

write_rds(models_stratum,
          str_c(folder_out, "models_stratum.rds"))
write_rds(model_combinations,
          str_c(folder_out, "model_combinations_stratum.rds"))
