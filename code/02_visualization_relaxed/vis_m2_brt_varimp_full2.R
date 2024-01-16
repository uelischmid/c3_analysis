## visualize boosted regression trees
## variable importance
## full models
## relaxed assessments
## 16.1.24, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)
library(dismo)

folder_in <- "data/processed/brt_models_v2_relaxed/"
folder_out <- "results/vis_models_v2_relaxed/brt/vis_brt_full_varimp2/"


# load data ---------------------------------------------------------------
models_full <- read_rds(str_c(folder_in, "models_full.rds"))
model_combinations <- read_rds(str_c(folder_in, "model_combinations_full.rds")) %>% 
  rowid_to_column()


# functions ---------------------------------------------------------------
plot_vi <- function(simt,
                    m_all  = models_full,
                    mc_all = model_combinations) {
  mc_used <- mc_all %>% 
    filter(simtype == simt)
  
  contributions <- vector(mode = "list", length = nrow(mc_used))
  
  for (i in 1:nrow(mc_used)) {
    vals <- mc_used[i, ]
    
    contributions[[i]] <- as_tibble(m_all[[pull(vals, rowid)]]$contributions) %>% 
      mutate(simtype       = pull(vals, simtype),
             nat_haz       = pull(vals, nat_haz),
             profile       = pull(vals, profile),
             resp_var_type = pull(vals, resp_var_type))
  }
  
  contributions <- bind_rows(contributions) %>% 
    mutate(var     = str_replace_all(var,
                                     c("stratum"       = "Zone",
                                       "q_reg2"        = "Qreg",
                                       "q_site2"       = "Qsite",
                                       "mgm_type"      = "Type",
                                       "mgm_interval"  = "Interval",
                                       "mgm_intensity" = "Intensity",
                                       "init"          = "Init")),
           var     = factor(var, levels = c("Zone", "Qsite", "Qreg", "Init",
                                            "Type", "Interval", "Intensity")),
           var     = fct_rev(var),
           nat_haz = factor(nat_haz, levels = c("A", "LED")),
           nat_haz = fct_rev(nat_haz),
           profile = factor(profile, levels = c("MP", "IP")),
           profile = fct_rev(profile))
  
  gg_p <- ggplot(contributions, aes(var, rel.inf, color = nat_haz, fill = profile)) +
    geom_col(position = "dodge",
             linewidth = 1) +
    coord_flip() +
    scale_color_discrete(type = c("#1b9e77", "#d95f02"),
                         guide = guide_legend(reverse = TRUE)) +
    scale_fill_discrete(type = c("#cccccc", "#525252"),
                        guide = guide_legend(reverse = TRUE)) +
    labs(x     = NULL,
         y     = "Variable importance (%)",
         color = "Natural\nhazard",
         fill  = "Profile") +
    facet_grid(cols = vars(resp_var_type),
               scales = "free_x") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  ggsave(filename = str_c(folder_out, "varimp_full_", simt, ".jpg"),
         plot     = gg_p,
         width = 16,
         height = 7,
         units = "cm",
         scale = 2)
}


# plot --------------------------------------------------------------------
plot_vi(simt = "ST")
plot_vi(simt = "LT")
