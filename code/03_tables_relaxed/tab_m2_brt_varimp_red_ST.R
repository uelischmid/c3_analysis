## table output of boosted regression trees
## variable importance
## reduced models ST
## relaxed assessments
## 2.11.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)
library(dismo)

folder_in <- "data/processed/brt_models_v2_relaxed/"
folder_out <- "results/tab_models_v2_relaxed/tab_brt_varimp/"

# load data ---------------------------------------------------------------
models <- read_rds(str_c(folder_in, "models_red_ST.rds"))
model_combinations <- read_rds(str_c(folder_in, "model_combinations_red_ST.rds")) %>% 
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
             stratum       = pull(vals, stratum),
             init          = pull(vals, init))
  }
  
  contributions <- bind_rows(contributions) %>% 
    mutate(var  = str_replace_all(var,
                                  c("q_site2"       = "Qsite",
                                    "mgm_type"      = "Type",
                                    "mgm_interval"  = "Interval",
                                    "mgm_intensity" = "Intensity")),
           init = str_replace_all(init,
                                  c("1" = "y",
                                    "2" = "s",
                                    "3" = "m")),
           init = factor(init, levels = c("y", "s", "m")))
  
  
  contributions <- contributions %>%
    mutate(pr_nh   = str_c(profile, "_", nat_haz),
           stratum = factor(stratum, levels = c("UM", "HM", "SA"))) %>% 
    dplyr::select(-c(nat_haz, profile)) %>% 
    pivot_wider(names_from = "pr_nh",
                values_from = "rel.inf") %>% 
    rowwise() %>% 
    mutate(vi_mean = mean(MP_A:IP_LED)) %>% 
    group_by(stratum, init) %>% 
    arrange(stratum, init, desc(vi_mean)) %>% 
    dplyr::select(simtype, resp_var_type, stratum, init, var, MP_A, MP_LED, IP_A, IP_LED, vi_mean)
  
  return(contributions)
}



# calculate ---------------------------------------------------------------
varimp_red_ST <- vector(mode = "list", length = 2)

varimp_red_ST[[1]] <- tab_vi(m_all  = models,
                             mc_all = model_combinations,
                             st     = "ST",
                             rv     = "abs")

varimp_red_ST[[2]] <- tab_vi(m_all  = models,
                             mc_all = model_combinations,
                             st     = "ST",
                             rv     = "diff")

varimp_red_ST <- bind_rows(varimp_red_ST)



write_delim(varimp_red_ST,
            str_c(folder_out, "varimp_red_ST.csv"),
            delim = ";")

