## visualize effects of betareg-models
## reduced models - diff
## all effects
## relaxed assessments
## differentiation between good and medium site quality
## 18.1.24, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(ggeffects)
library(cowplot)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v2_relaxed/"
folder_out <- "results/vis_models_v2_relaxed/betareg/vis_betareg_red_effects_combined_4b/"


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
  effs_type_g <- get_eff(mod  = m_all,
                         mc   = mc_all %>% 
                           filter(stratum == stra) %>% 
                           filter(init == st_init) %>% 
                           filter(resp_var_type == "diff"),
                         term = c("mgm_type",
                                  "q_site2 [good]")) %>% 
    mutate(profile  = factor(profile, levels = c("MP", "IP")))
  
  effs_type_m <- get_eff(mod  = m_all,
                         mc   = mc_all %>% 
                           filter(stratum == stra) %>% 
                           filter(init == st_init) %>% 
                           filter(resp_var_type == "diff"),
                         term = c("mgm_type",
                                  "q_site2 [medium]")) %>% 
    mutate(profile  = factor(profile, levels = c("MP", "IP")))
  
  
  effs_intint_g <- get_eff(mod  = m_all,
                           mc   = mc_all %>% 
                             filter(stratum == stra) %>% 
                             filter(init == st_init) %>% 
                             filter(resp_var_type == "diff"),
                           term = c("mgm_interval [10:40]",
                                    "mgm_intensity [10:40]",
                                    "q_site2 [good]")) %>% 
    rename(interval  = x,
           intensity = group) %>% 
    mutate(intensity = as.character(intensity),
           intensity = as.numeric(intensity),
           profile  = factor(profile, levels = c("MP", "IP")))
  
  effs_intint_m <- get_eff(mod  = m_all,
                           mc   = mc_all %>% 
                             filter(stratum == stra) %>% 
                             filter(init == st_init) %>% 
                             filter(resp_var_type == "diff"),
                           term = c("mgm_interval [10:40]",
                                    "mgm_intensity [10:40]",
                                    "q_site2 [medium]")) %>% 
    rename(interval  = x,
           intensity = group) %>% 
    mutate(intensity = as.character(intensity),
           intensity = as.numeric(intensity),
           profile  = factor(profile, levels = c("MP", "IP")))
  
  
  ylims_type <- c(min(c(effs_type_g$conf.low, effs_type_m$conf.low)),
                  max(c(effs_type_g$conf.high, effs_type_m$conf.high)))
  if (ylims_type[1] > 0) ylims_type[1] <- 0
  if (ylims_type[2] < 0) ylims_type[2] <- 0
  
  ylims_intint <- c(min(c(effs_intint_g$predicted, effs_intint_m$predicted)),
                    max(c(effs_intint_g$predicted, effs_intint_m$predicted)))
  if (ylims_intint[1] > 0) ylims_intint[1] <- 0
  if (ylims_intint[2] < 0) ylims_intint[2] <- 0
  
  # plot titles
  gg_title_g <- ggplot() +
    annotate(geom     = "text",
             x        = 0,
             y        = 0.5,
             label    = expression(Q[site]*": good"),
             angle    = 0,
             fontface = "bold") +
    theme_void()
  
  gg_title_m <- ggplot() +
    annotate(geom     = "text",
             x        = 0,
             y        = 0.5,
             label    = expression(Q[site]*": medium"),
             angle    = 0,
             fontface = "bold") +
    theme_void()
  
  # plot type
  gg_type_g_temp <- ggplot(effs_type_g, aes(x, predicted, color = nat_haz)) +
    geom_hline(yintercept = 0,
               color      = "darkgrey") +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    scale_y_continuous(limits = ylims_type) +
    labs(x     = "Type",
         y     = "\u0394 PQ",
         color = "Natural\nhazard") +
    facet_grid(cols = vars(profile)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.text.x      = element_text(angle = 45, hjust = 1))
  
  gg_type_g <- gg_type_g_temp +
    theme(legend.position = "none")
  
  legend_nathaz <- get_legend(gg_type_g_temp)
  
  
  gg_type_m <- ggplot(effs_type_m, aes(x, predicted, color = nat_haz)) +
    geom_hline(yintercept = 0,
               color      = "darkgrey") +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    scale_y_continuous(limits = ylims_type) +
    labs(x     = "Type",
         y     = "\u0394 PQ",
         color = "Natural\nhazard") +
    facet_grid(cols = vars(profile)) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position  = "none",
          axis.text.x      = element_text(angle = 45, hjust = 1))
  
  # plot heatmap intint
  gg_intint_g_temp <- ggplot(effs_intint_g, aes(interval, intensity, fill = predicted)) +
    geom_tile() +
    scale_fill_gradient2(low      = "#ff0000",
                         mid      = "#ffffff",
                         high     = "#004a8d",
                         midpoint = 0,
                         limits   = ylims_intint) +
    labs(x     = "Interval (y)",
         y     = "Intensity (%)",
         fill  = "\u0394 PQ") +
    facet_grid(rows = vars(nat_haz),
               cols = vars(profile)) +
    coord_equal() +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  gg_intint_g <- gg_intint_g_temp +
    theme(legend.position = "none")
  
  legend_deltapq <- get_legend(gg_intint_g_temp)
  
  gg_intint_m <- ggplot(effs_intint_m, aes(interval, intensity, fill = predicted)) +
    geom_tile() +
    scale_fill_gradient2(low      = "#ff0000",
                         mid      = "#ffffff",
                         high     = "#004a8d",
                         midpoint = 0,
                         limits   = ylims_intint) +
    labs(x     = "Interval (y)",
         y     = "Intensity (%)",
         fill  = "\u0394 PQ") +
    facet_grid(rows = vars(nat_haz),
               cols = vars(profile)) +
    coord_equal() +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # plot all
  plots <- list()
  plots[[1]] <- gg_title_g
  plots[[2]] <- gg_title_m
  plots[[4]] <- gg_type_g
  plots[[5]] <- gg_type_m
  plots[[6]] <- legend_nathaz
  plots[[7]] <- gg_intint_g
  plots[[8]] <- gg_intint_m
  plots[[9]] <- legend_deltapq
  
  gg_out <- plot_grid(plotlist = plots,
                      ncol = 3,
                      rel_heights = c(0.2, 1, 1.5),
                      rel_widths = c(1, 1, 0.2),
                      labels = c("", "", "", "a", "b", "", "c", "d", ""),
                      align = "v",
                      axis = c("lr"))
  
  ggsave(filename = str_c(f_out, "eff_comb_diff_",
                          stra, "_", st_init,
                          "_qsite.jpg"),
         plot = gg_out,
         width = 15,
         height = 7.5,
         units = "cm",
         scale = 1.8)
}

# plot --------------------------------------------------------------------
plot_vars <- expand_grid(stratum    = c("UM", "HM", "SA"),
                         stand_init = c("1", "2", "3", "LT"))

for (i in 1:nrow(plot_vars)) {
  vars <- plot_vars[i, ]
  plot_effs(stra    = pull(vars, stratum),
            st_init = pull(vars, stand_init))
}
