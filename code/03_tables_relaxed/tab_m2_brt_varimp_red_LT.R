## table output of boosted regression trees
## variable importance
## reduced models LT
## relaxed assessments
## 2.11.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)
library(dismo)

folder_in <- "data/processed/brt_models_v2_relaxed/"
folder_out <- "results/tab_models_v2_relaxed/tab_brt_varimp/"

# load data ---------------------------------------------------------------
models <- read_rds(str_c(folder_in, "models_red_LT.rds"))
model_combinations <- read_rds(str_c(folder_in, "model_combinations_red_LT.rds")) %>% 
  rowid_to_column()


# functions ---------------------------------------------------------------
tab_vi <- function(m_all, mc_all,
                   st, rv) {
  mc_used <- mc_all %>% 
    filter(simtype == st) %>% 
    filter(resp_var_type == rv)
  
  contributions <- vector(mode = "list", length = nrow(mc_used))
  
  for (i in 1:nrow(mc_used)) {
    vals <- mc_used[i, ]
    
    contributions[[i]] <- as_tibble(m_all[[pull(vals, rowid)]]$contributions) %>% 
      mutate(simtype       = pull(vals, simtype),
             nat_haz       = pull(vals, nat_haz),
             profile       = pull(vals, profile),
             resp_var_type = pull(vals, resp_var_type),
             stratum       = pull(vals, stratum))
  }
  
  contributions <- bind_rows(contributions) %>% 
    mutate(var   = str_replace_all(var,
                                   c("q_site2"       = "Qsite",
                                     "mgm_type"      = "Type",
                                     "mgm_interval"  = "Interval",
                                     "mgm_intensity" = "Intensity")))
  
  
  contributions <- contributions %>%
    mutate(pr_nh   = str_c(profile, "_", nat_haz),
           stratum = factor(stratum, levels = c("UM", "HM", "SA"))) %>% 
    dplyr::select(-c(nat_haz, profile)) %>% 
    pivot_wider(names_from = "pr_nh",
                values_from = "rel.inf") %>% 
    rowwise() %>% 
    mutate(vi_mean = mean(MP_A:IP_LED)) %>% 
    group_by(stratum) %>% 
    arrange(stratum, desc(vi_mean)) %>% 
    dplyr::select(simtype, resp_var_type, stratum, var, MP_A, MP_LED, IP_A, IP_LED, vi_mean)
  
  return(contributions)
}



# calculate ---------------------------------------------------------------
varimp_red_LT <- vector(mode = "list", length = 2)

varimp_red_LT[[1]] <- tab_vi(m_all  = models,
                             mc_all = model_combinations,
                             st     = "LT",
                             rv     = "abs")

varimp_red_LT[[2]] <- tab_vi(m_all  = models,
                             mc_all = model_combinations,
                             st     = "LT",
                             rv     = "diff")

varimp_red_LT <- bind_rows(varimp_red_LT)



write_delim(varimp_red_LT,
            str_c(folder_out, "varimp_red_LT.csv"),
            delim = ";")

