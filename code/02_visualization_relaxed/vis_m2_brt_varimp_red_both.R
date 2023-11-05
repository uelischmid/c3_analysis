## visualize boosted regression trees
## variable importance
## reduced models
## relaxed assessments
## 5.11.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)
library(dismo)

folder_in <- "data/processed/brt_models_v2_relaxed/"
folder_out <- "results/vis_models_v2_relaxed/brt/vis_brt_red_varimp/"

# load data ---------------------------------------------------------------
ST_models <- read_rds(str_c(folder_in, "models_red_ST.rds"))
ST_model_combinations <- read_rds(str_c(folder_in, "model_combinations_red_ST.rds")) %>% 
  rowid_to_column()

LT_models <- read_rds(str_c(folder_in, "models_red_LT.rds"))
LT_model_combinations <- read_rds(str_c(folder_in, "model_combinations_red_LT.rds")) %>% 
  rowid_to_column() %>% 
  mutate(init = "LT")

# functions ---------------------------------------------------------------
get_varimp <- function(rv, m_all, mc_all) {
  mc_used <- mc_all %>% 
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
    mutate(var   = str_replace_all(var,
                                   c("q_site2"       = "Qsite",
                                     "mgm_type"      = "Type",
                                     "mgm_interval"  = "Interval",
                                     "mgm_intensity" = "Intensity")),
           init2 = case_when(init == "1" ~ "ST: young",
                             init == "2" ~ "ST: structured",
                             init == "3" ~ "ST: mature",
                             TRUE        ~ "LT"))
  
  return(contributions)
}


# get variable importances ------------------------------------------------
# diff
varimp_ST_diff <- get_varimp(rv     = "diff",
                             m_all  = ST_models,
                             mc_all = ST_model_combinations)

varimp_LT_diff <- get_varimp(rv     = "diff",
                             m_all  = LT_models,
                             mc_all = LT_model_combinations)

varimp_both_diff <- bind_rows(varimp_ST_diff, varimp_LT_diff) %>% 
  mutate(var     = factor(var, levels = c("Type", "Interval", "Intensity", "Qsite")),
         var     = fct_rev(var),
         nat_haz = factor(nat_haz, levels = c("A", "LED")),
         nat_haz = fct_rev(nat_haz),
         profile = factor(profile, levels = c("MP", "IP")),
         profile = fct_rev(profile),
         stratum = factor(stratum, levels = c("UM", "HM", "SA")),
         init2   = factor(init2, levels = c("ST: young", "ST: structured", "ST: mature", "LT")))



# plot --------------------------------------------------------------------
ggplot(varimp_both_diff, aes(var, rel.inf, color = nat_haz, fill = profile)) +
  geom_col(position = "dodge",
           linewidth = 0.5) +
  coord_flip() +
  scale_color_discrete(type = c("#1b9e77", "#d95f02"),
                       guide = guide_legend(reverse = TRUE)) +
  scale_fill_discrete(type = c("#cccccc", "#525252"),
                      guide = guide_legend(reverse = TRUE)) +
  facet_grid(rows     = vars(init2),
             cols     = vars(stratum)) +
  labs(x     = NULL,
       y     = "Variable importance (%)",
       color = "Natural\nhazard",
       fill  = "Profile") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(filename = str_c(folder_out, "both_diff.png"),
       width    = 20,
       height   = 20,
       units    = "cm")


