## visualize boosted regression trees
## variable contributions
## 23.2.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)
library(dismo)
library(cowplot)

folder_in <- "data/processed/brt_models_v1/"
folder_out <- "results/vis_models_v1/brt/vis_brt_varcont_stratum/"


# load data ---------------------------------------------------------------
models_full <- read_rds(str_c(folder_in, "models_full.rds"))
models_stratum <- read_rds(str_c(folder_in, "models_stratum.rds"))
model_combinations_full <- read_rds(str_c(folder_in, "model_combinations_full.rds")) %>% 
  rowid_to_column()
model_combinations_stratum <- read_rds(str_c(folder_in, "model_combinations_stratum.rds")) %>% 
  rowid_to_column()



# plot --------------------------------------------------------------------
for(a in 1:nrow(model_combinations_full)) {
  mc_full <- model_combinations_full[a, ]
  mc_stratum <- model_combinations_stratum %>% 
    filter(simtype == pull(mc_full, simtype)) %>% 
    filter(nat_haz == pull(mc_full, nat_haz)) %>% 
    filter(profile == pull(mc_full, profile)) %>% 
    filter(resp_var_type == pull(mc_full, resp_var_type))
  
  cont_full <- as_tibble(models_full[[pull(mc_full, rowid)]]$contributions) %>% 
    mutate(stratum = "combined")
  
  cont_str <- vector(mode = "list", length = nrow(mc_stratum))
  for (i in 1:nrow(mc_stratum)) {
    cont_str[[i]] <- as_tibble(models_stratum[[pull(mc_stratum, rowid)[i]]]$contributions) %>% 
      mutate(stratum = pull(mc_stratum, stratum)[i])
  }
  cont_str <- bind_rows(cont_str)
  
  cont_all <- bind_rows(cont_full, cont_str) %>% 
    mutate(var = str_replace_all(var,
                                 c("stratum"       = "Stratum",
                                   "q_reg2"        = "Qreg",
                                   "q_site2"       = "Qsite",
                                   "mgm_type"      = "Type",
                                   "mgm_interval"  = "Interval",
                                   "mgm_intensity" = "Intensity",
                                   "init"          = "Init")))
  var_order <- cont_all %>% 
    filter(stratum == "combined") %>% 
    arrange(desc(rel.inf)) %>% 
    pull(var)
  cont_all <- cont_all %>% 
    mutate(var     = factor(var, levels = var_order),
           var     = fct_rev(var),
           stratum = factor(stratum, levels = c("combined", "UM", "HM", "SA")))
  
  cont_mgm <- cont_all %>% 
    mutate(mgm_var = case_when(var %in% c("Type", "Interval", "Intensity") ~ "mgm",
                               TRUE                                        ~ "other")) %>% 
    group_by(stratum, mgm_var) %>% 
    summarise(rel.inf.sum = sum(rel.inf),
              .groups     = "drop") %>% 
    mutate(stratum = fct_rev(stratum))
  
  title_str <- str_c(pull(mc_full, simtype),
                     pull(mc_full, nat_haz),
                     pull(mc_full, profile),
                     pull(mc_full, resp_var_type),
                     sep = " ")
  # all variables
  p1 <- ggplot(cont_all, aes(var, rel.inf)) +
    geom_col() +
    coord_flip() +
    facet_wrap(~stratum,
               nrow = 1,
               scales = "free_x") +
    labs(title = title_str,
         x     = NULL,
         y     = "Variable importance (%)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  ggsave(filename = str_c(folder_out,
                          str_replace_all(title_str, " ", "_"),
                          ".png"),
         plot     = p1)
  
  # mgm
  p2 <- cont_mgm %>% 
    filter(mgm_var == "mgm") %>% 
    ggplot(aes(stratum, rel.inf.sum)) +
    geom_col() +
    scale_y_continuous(limits = c(0, 100)) +
    coord_flip() +
    labs(title = title_str,
         x     = NULL,
         y     = "Importance of management (%)") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  ggsave(filename = str_c(folder_out,
                          str_replace_all(title_str, " ", "_"),
                          "_mgm.png"),
         plot     = p2)
}



