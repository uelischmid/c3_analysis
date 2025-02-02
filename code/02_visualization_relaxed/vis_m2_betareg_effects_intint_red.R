## visualize effects of betareg-models
## reduced models - diff
## intensity & interval
## relaxed assessments
## 9.11.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(ggeffects)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v2_relaxed/"
folder_out <- "results/vis_models_v2_relaxed/betareg/vis_betareg_red_effects_intint/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in_d, "analysis_data_transf_relaxed.rds")) %>% 
  filter(q_reg2 == "normal")


# subset ST data and save as new objects
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
    filter(init    == pull(vals, init)) %>% 
    mutate(sha_i_MP_met_abs_t = (sha_i_MP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.),
           sha_i_IP_met_abs_t = (sha_i_IP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.))
  
  assign(object_name, object_content)
}

rm(object_name, object_content)

# subset LT data and save as new objects
data_comb <- expand_grid(simtype = "LT",
                         nat_haz = c("A", "LED"),
                         stratum = c("UM", "HM", "SA"))
for (i in 1:nrow(data_comb)) {
  vals <- data_comb[i,]
  object_name <- str_c("data_",
                       pull(vals, simtype), "_",
                       pull(vals, nat_haz), "_",
                       pull(vals, stratum))
  object_content <- data_analysis %>% 
    filter(simtype == pull(vals, simtype)) %>% 
    filter(nat_haz == pull(vals, nat_haz)) %>% 
    filter(stratum == pull(vals, stratum)) %>% 
    mutate(sha_i_MP_met_abs_t = (sha_i_MP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.),
           sha_i_IP_met_abs_t = (sha_i_IP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.))
  
  assign(object_name, object_content)
}

rm(object_name, object_content)


# load models -------------------------------------------------------------
models <- c(read_rds(str_c(folder_in_m, "models_red_ST.rds")),
            read_rds(str_c(folder_in_m, "models_red_LT.rds")))
model_combinations <- bind_rows(read_rds(str_c(folder_in_m, "model_combinations_red_ST.rds")),
                                mc_LT <- read_rds(str_c(folder_in_m, "model_combinations_red_LT.rds")) %>% 
                                  mutate(init = "LT")) %>% 
  rowid_to_column()



# functions ---------------------------------------------------------------
get_eff <- function(mod, mc, term) {
  effs <- vector(mode = "list", length = nrow(mc))
  
  for (i in 1:nrow(mc)) {
    vals <- mc[i, ]
    eff <- ggeffect(model = mod[[pull(vals, rowid)]],
                    terms = term) %>% 
      as_tibble() %>% 
      mutate(nat_haz = pull(vals, nat_haz),
             profile = pull(vals, profile),
             init    = pull(vals, init))
    
    # backtransform for diff-models
    if (pull(vals, resp_var_type) == "diff") {
      eff <- eff %>% 
        mutate(predicted = 2 * predicted - 1,
               std.error = 2 * std.error - 1,
               conf.low  = 2 * conf.low - 1,
               conf.high = 2 * conf.high - 1)
    }
    
    effs[[i]] <- eff
  }
  effs <- bind_rows(effs)
  return(effs)
}

plot_heatmap_intint <- function(stra, pro,
                                m_all  = models,
                                mc_all = model_combinations,
                                f_out  = folder_out) {
  effs <- get_eff(mod  = m_all,
                  mc   = mc_all %>% 
                    filter(stratum == stra) %>% 
                    filter(profile == pro) %>% 
                    filter(resp_var_type == "diff"),
                  term = c("mgm_interval [10:40]",
                           "mgm_intensity [10:40]")) %>% 
    rename(interval  = x,
           intensity = group) %>% 
    mutate(intensity = as.character(intensity),
           intensity = as.numeric(intensity),
           init      = factor(init),
           init      = fct_recode(init,
                                  "ST: young"      = "1",
                                  "ST: structured" = "2",
                                  "ST: mature"     = "3"))
  
  gg <- ggplot(effs, aes(interval, intensity, fill = predicted)) +
    geom_tile() +
    scale_fill_gradient2(low      = "#ff0000",
                         mid      = "#ffffff",
                         high     = "#004a8d",
                         midpoint = 0) +
    labs(x     = "Interval (y)",
         y     = "Intensity (%)",
         fill  = "\u0394 PQ") +
    facet_grid(rows = vars(nat_haz),
               cols = vars(init)) +
    coord_equal() +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  ggsave(filename = str_c(f_out, "heatmap_intint_",
                          stra, "_",
                          pro, "_diff.jpg"),
         plot     = gg,
         width    = 16,
         height   = 8,
         units    = "cm",
         scale    = 1)
}

# plot --------------------------------------------------------------------
plot_vars <- expand_grid(stratum = c("UM", "HM", "SA"),
                         profile = c("MP", "IP"))


for (i in 1:nrow(plot_vars)) {
  vars <- plot_vars[i, ]
  plot_heatmap_intint(stra = pull(vars, stratum),
                      pro  = pull(vars, profile))
}





