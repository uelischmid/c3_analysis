## visualize boosted regression trees
## variable importance
## reduced models LT
## relaxed assessments
## 12.5.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)
library(dismo)

folder_in <- "data/processed/brt_models_v2_relaxed/"
folder_out <- "results/vis_models_v2_relaxed/brt/vis_brt_red_varimp/"

# load data ---------------------------------------------------------------
models <- read_rds(str_c(folder_in, "models_red_ST.rds"))
model_combinations <- read_rds(str_c(folder_in, "model_combinations_red_ST.rds")) %>% 
  rowid_to_column()



# functions ---------------------------------------------------------------
plot_vi <- function(st, rv,
                    m_all  = models,
                    mc_all = model_combinations,
                    f_out  = folder_out) {
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
    mutate(var   = str_replace_all(var,
                                   c("q_site2"       = "Qsite",
                                     "mgm_type"      = "Type",
                                     "mgm_interval"  = "Interval",
                                     "mgm_intensity" = "Intensity")),
           init2 = case_when(init == "1" ~ "young",
                             init == "2" ~ "structured",
                             TRUE        ~ "mature"))
  
  # var_order <- contributions %>% 
  #   group_by(var) %>% 
  #   summarise(m_inf = mean(rel.inf)) %>% 
  #   arrange(desc(m_inf)) %>% 
  #   pull(var)
  
  contributions <- contributions %>% 
    mutate(# var     = factor(var, levels = var_order),
           var     = factor(var, levels = c("Type", "Interval", "Intensity", "Qsite")),
           var     = fct_rev(var),
           nat_haz = factor(nat_haz, levels = c("A", "LED")),
           nat_haz = fct_rev(nat_haz),
           profile = factor(profile, levels = c("MP", "IP")),
           profile = fct_rev(profile),
           stratum = factor(stratum, levels = c("UM", "HM", "SA")),
           init2   = factor(init2, levels = c("young", "structured", "mature")))
  
  title_str <- str_c(st, " - ", rv)
  
  gg_p <- ggplot(contributions, aes(var, rel.inf, color = nat_haz, fill = profile)) +
    geom_col(position = "dodge",
             linewidth = 0.5) +
    coord_flip() +
    scale_color_discrete(type = c("#1b9e77", "#d95f02"),
                         guide = guide_legend(reverse = TRUE)) +
    scale_fill_discrete(type = c("#cccccc", "#525252"),
                        guide = guide_legend(reverse = TRUE)) +
    facet_grid(rows     = vars(init2),
               cols     = vars(stratum)) +
    labs(title = title_str,
         x     = NULL,
         y     = "Variable importance (%)",
         color = "Natural\nhazard",
         fill  = "Profile") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  ggsave(filename = str_c(f_out, st, "_", rv, ".png"),
         plot     = gg_p,
         width    = 20,
         height   = 15,
         units    = "cm")
}


# plot --------------------------------------------------------------------
# ST abs
plot_vi(st = "ST",
        rv = "abs")

# ST diff
plot_vi(st = "ST",
        rv = "diff")

