## visualize effects of betareg-models
## reduced models - diff
## all effects
## relaxed assessments
## 15.1.24, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(ggeffects)
library(cowplot)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v2_relaxed/"
folder_out <- "results/vis_models_v2_relaxed/betareg/vis_betareg_red_effects_combined_3/"


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

plot_effs <- function(stra, st_init,
                      m_all  = models,
                      mc_all = model_combinations,
                      f_out  = folder_out) {
  # effects
  effs_qsite <- get_eff(mod  = m_all,
                        mc   = mc_all %>% 
                          filter(stratum == stra) %>% 
                          filter(init == st_init) %>% 
                          filter(resp_var_type == "diff"),
                        term = c("q_site2")) %>% 
    mutate(profile  = factor(profile, levels = c("MP", "IP")))
  
  effs_type <- get_eff(mod  = m_all,
                       mc   = mc_all %>% 
                         filter(stratum == stra) %>% 
                         filter(init == st_init) %>% 
                         filter(resp_var_type == "diff"),
                       term = c("mgm_type")) %>% 
    mutate(profile  = factor(profile, levels = c("MP", "IP")))
  
  
  effs_intint <- get_eff(mod  = m_all,
                         mc   = mc_all %>% 
                           filter(stratum == stra) %>% 
                           filter(init == st_init) %>% 
                           filter(resp_var_type == "diff"),
                         term = c("mgm_interval [10:40]",
                                  "mgm_intensity [10:40]")) %>% 
    rename(interval  = x,
           intensity = group) %>% 
    mutate(intensity = as.character(intensity),
           intensity = as.numeric(intensity),
           profile  = factor(profile, levels = c("MP", "IP")))
  
  
  ylims_upper <- c(min(c(effs_qsite$conf.low, effs_type$conf.low)),
                   max(c(effs_qsite$conf.high, effs_type$conf.high)))
  if (ylims_upper[1] > 0) ylims_upper[1] <- 0
  if (ylims_upper[2] < 0) ylims_upper[2] <- 0
  
  # plot qsite
  gg_qsite_temp <- ggplot(effs_qsite, aes(x, predicted, color = nat_haz)) +
    geom_hline(yintercept = 0,
               color      = "darkgrey") +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    scale_y_continuous(limits = ylims_upper) +
    labs(x = expression(Q[site]),
         y = "\u0394 PQ",
         color = "Natural\nhazard") +
    facet_grid(cols = vars(profile)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  gg_qsite <- gg_qsite_temp +
    theme(legend.position = "none")
  
  legend_nathaz <- get_legend(gg_qsite_temp)
  
  # plot type
  gg_type <- ggplot(effs_type, aes(x, predicted, color = nat_haz)) +
    geom_hline(yintercept = 0,
               color      = "darkgrey") +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    scale_y_continuous(limits = ylims_upper) +
    labs(x = "Type",
         y = "\u0394 PQ") +
    facet_grid(cols = vars(profile)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank()) +
    theme(legend.position = "none")
  
  # plot heatmap intint
  gg_intint_temp <- ggplot(effs_intint, aes(interval, intensity, fill = predicted)) +
    geom_tile() +
    scale_fill_gradient2(low      = "#ff0000",
                         mid      = "#ffffff",
                         high     = "#004a8d",
                         midpoint = 0) +
    labs(x     = "Interval (y)",
         y     = "Intensity (%)",
         fill  = "\u0394 PQ") +
    facet_grid(rows = vars(nat_haz),
               cols = vars(profile)) +
    coord_equal() +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  gg_intint <- gg_intint_temp +
    theme(legend.position = "none")
  
  legend_deltapq <- get_legend(gg_intint_temp)
  
  # plot all
  plots <- list()
  plots[[1]] <- gg_qsite
  plots[[2]] <- legend_nathaz
  plots[[3]] <- gg_type
  plots[[5]] <- gg_intint
  plots[[6]] <- legend_deltapq
  
  gg_out <- plot_grid(plotlist = plots,
            ncol = 2,
            rel_heights = c(1, 1, 2),
            rel_widths = c(1, 0.2),
            labels = c("a", "", "b", "", "c", ""),
            align = "v",
            axis = c("lr"))
  
  ggsave(filename = str_c(f_out, "eff_comb_diff_",
                          stra, "_", st_init, ".jpg"),
         plot = gg_out,
         width = 8,
         height = 10,
         units = "cm",
         scale = 2)
}

# plot --------------------------------------------------------------------
plot_vars <- expand_grid(stratum    = c("UM", "HM", "SA"),
                         stand_init = c("1", "2", "3", "LT"))

for (i in 1:nrow(plot_vars)) {
  vars <- plot_vars[i, ]
  plot_effs(stra    = pull(vars, stratum),
            st_init = pull(vars, stand_init))
}
