## visualize effects of betareg-models
## full models
## relaxed assessments
## effects of and with qreg
## 9.1.24, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(ggeffects)
library(cowplot)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v2_relaxed/"
folder_out <- "results/vis_models_v2_relaxed/betareg/vis_betareg_full_effects_combined_qreg/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in_d, "analysis_data_transf_relaxed.rds"))

data_LT_A <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "A") %>% 
  mutate(sha_i_MP_met_abs_t = (sha_i_MP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.),
         sha_i_IP_met_abs_t = (sha_i_IP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.))

data_LT_LED <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "LED") %>% 
  mutate(sha_i_MP_met_abs_t = (sha_i_MP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.),
         sha_i_IP_met_abs_t = (sha_i_IP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.))

data_ST_A <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "A") %>% 
  mutate(init = fct_drop(init)) %>% 
  mutate(sha_i_MP_met_abs_t = (sha_i_MP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.),
         sha_i_IP_met_abs_t = (sha_i_IP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.))

data_ST_LED <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "LED") %>% 
  mutate(init = fct_drop(init)) %>% 
  mutate(sha_i_MP_met_abs_t = (sha_i_MP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.),
         sha_i_IP_met_abs_t = (sha_i_IP_met_abs * (nrow(.)- 1) + 0.5) / nrow(.))


# load models -------------------------------------------------------------
models <- read_rds(str_c(folder_in_m, "models_full.rds"))
model_combinations <- read_rds(str_c(folder_in_m, "model_combinations_full.rds")) %>% 
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
             rvt     = pull(vals, resp_var_type),
             simtype = pull(vals, simtype))
    
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


plot_profilename <- function(name) {
  gg_p <- ggplot() +
    annotate(geom     = "text",
             x        = 0,
             y        = 0.5,
             label    = name,
             angle    = 0,
             fontface = "bold") +
    theme_void()
  return(gg_p)
}

plot_predictorname <- function(name) {
  gg_p <- ggplot() +
    annotate(geom     = "text",
             x        = 0,
             y        = 0.5,
             label    = name,
             angle    = 90,
             fontface = "bold") +
    theme_void()
  return(gg_p)
}


plot_effects_comb_mgm_qreg <- function(pr, rvt,
                                       m_all  = models,
                                       mc_all = model_combinations,
                                       f_out  = folder_out) {
  # select model combinations
  mc_used <- mc_all %>% 
    filter(profile == pr) %>% 
    filter(resp_var_type == rvt)
  
  # calculate effects
  effs_mgm_type <- get_eff(mod  = m_all,
                           mc   = mc_used,
                           term = c("mgm_type",
                                    "q_reg2"))
  effs_mgm_intint <- get_eff(mod  = m_all,
                             mc   = mc_used,
                             term = c("mgm_interval [10:40]",
                                      "mgm_intensity [10, 20, 30, 40]",
                                      "q_reg2"))
  
  # prepare plots
  theme_custom <- theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  effs_all <- bind_rows(effs_mgm_type,
                        effs_mgm_intint %>% 
                          mutate(x = factor(x)))
  
  y_lims <- c(floor(min(effs_all$conf.low) * 20) / 20,
             ceiling(max(effs_all$conf.high) * 20) / 20)
  
  if (y_lims[2] - y_lims[1] < 0.1) {
    y_breaks <- seq(y_lims[1], y_lims[2], 0.02)
  } else if (y_lims[2] - y_lims[1] < 0.3) {
    y_breaks <- seq(y_lims[1], y_lims[2], 0.05)
  } else {
    y_breaks <- seq(y_lims[1], y_lims[2], 0.1)
  }
  
  if (rvt == "abs") {
    label_y <- "Profile met (%)"
    scale_y_custom_ST <- scale_y_continuous(limits = y_lims,
                                            breaks = y_breaks,
                                            labels = scales::percent)
    scale_y_custom_LT <- scale_y_continuous(limits = y_lims,
                                            breaks = y_breaks,
                                            labels = NULL)
    
  } else {
    label_y <- "\u0394 PQ"
    scale_y_custom_ST <- scale_y_continuous(limits = y_lims,
                                            breaks = y_breaks)
    scale_y_custom_LT <- scale_y_continuous(limits = y_lims,
                                            breaks = y_breaks,
                                            labels = NULL)
  }
  
  # name plots
  p_title_ST <- plot_profilename("Short-term simulations")
  p_title_LT <- plot_profilename("Long-term simulations")
  p_title_type <- plot_predictorname("Type")
  p_title_intint <- plot_predictorname("Interval * Intensity")

  # plot mgm type
  # ST
  gg_type_ST <- effs_mgm_type %>% 
    rename(Qreg = group) %>% 
    filter(simtype == "ST") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_hline(yintercept = 0,
               color      = "darkgrey") +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x     = "Type",
         y     = label_y,
         color = "Natural\nhazard") +
    facet_grid(rows = vars(Qreg),
               labeller = "label_both") +
    scale_y_custom_ST +
    theme_custom
    
  # LT
  gg_type_LT <- effs_mgm_type %>% 
    rename(Qreg = group) %>% 
    filter(simtype == "LT") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_hline(yintercept = 0,
               color      = "darkgrey") +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x     = "Type",
         y     = "",
         color = "Natural\nhazard") +
    facet_grid(rows = vars(Qreg),
               labeller = "label_both") +
    scale_y_custom_LT +
    theme_custom
  
  
  # plot intensity * interval interaction (x = interval)
  # ST
  gg_intint_ST <- effs_mgm_intint %>% 
    rename(Intensity = group) %>% 
    filter(simtype == "ST") %>% 
    mutate(nhfa = str_c(nat_haz, facet)) %>% 
    ggplot(aes(x, predicted, color = nat_haz, lty = facet)) +
    geom_hline(yintercept = 0, color = "darkgrey") +
    geom_ribbon(aes(x, predicted,
                    ymin  = conf.low,
                    ymax  = conf.high,
                    group = nhfa,
                    fill  = nat_haz),
                alpha       = 0.2,
                color       = NA,
                inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(~Intensity,
               labeller = "label_both") +
    labs(x     = "Interval (y)",
         y     = label_y,
         color = "Natural\nhazard",
         fill  = "Natural\nhazard",
         lty   = expression(Q[reg])) +
    scale_y_custom_ST +
    theme_custom
    
  # LT
  gg_intint_LT <- effs_mgm_intint %>% 
    rename(Intensity = group) %>% 
    filter(simtype == "LT") %>% 
    mutate(nhfa = str_c(nat_haz, facet)) %>% 
    ggplot(aes(x, predicted, color = nat_haz, lty = facet)) +
    geom_hline(yintercept = 0, color = "darkgrey") +
    geom_ribbon(aes(x, predicted,
                    ymin  = conf.low,
                    ymax  = conf.high,
                    group = nhfa,
                    fill  = nat_haz),
                alpha       = 0.2,
                color       = NA,
                inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(~Intensity,
               labeller = "label_both") +
    labs(x     = "Interval (y)",
         y     = "",
         color = "Natural\nhazard",
         fill  = "Natural\nhazard",
         lty   = expression(Q[reg])) +
    scale_y_custom_LT +
    theme_custom
  
  
  # plot together
  plots <- vector(mode = "list", length = 12)
  
  plots[[2]] <- p_title_ST
  plots[[3]] <- p_title_LT
  
  plots[[5]] <- p_title_type
  plots[[6]] <- gg_type_ST +
    theme(legend.position = "none")
  plots[[7]] <- gg_type_LT +
    theme(legend.position = "none")
  plots[[8]] <- get_legend(gg_type_ST)
  
  plots[[9]] <- p_title_intint
  plots[[10]] <- gg_intint_ST +
    theme(legend.position = "none")
  plots[[11]] <- gg_intint_LT +
    theme(legend.position = "none")
  plots[[12]] <- get_legend(gg_intint_ST)
  
  
  gg_out <- plot_grid(plotlist    = plots,
                      ncol        = 4,
                      rel_heights = c(0.15, 1, 1),
                      rel_widths  = c(0.15, 1, 0.9, 0.4),
                      align       = "h",
                      axis        = "lr")
  
  ggsave(filename = str_c(f_out, "betareg_full_mgm_",
                          pr, "_",
                          rvt, ".jpg"),
         plot     = gg_out,
         width    = 16,
         height   = 12,
         units    = "cm",
         scale    = 1.5)
}


# plot qreg ---------------------------------------------------------------
# calculate effects
effs_qual_reg <- get_eff(mod  = models,
                         mc   = model_combinations,
                         term = "q_reg2") %>% 
  mutate(simtype    = factor(simtype, levels = c("ST", "LT")),
         profile    = factor(profile, levels = c("MP", "IP")))

# prepare plots
theme_custom <- theme_bw() +
  theme(panel.grid.minor = element_blank())

# plot regeneration quality
# plot abs
p_qreg_abs <- effs_qual_reg %>% 
  filter(rvt == "abs") %>% 
  ggplot(aes(x, predicted, color = nat_haz)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width    = 0.25,
                position = position_dodge(width = 0.2)) +
  labs(x     = expression(Q[reg]),
       y     = "Profile met (%)",
       color = "Natural\nhazard") +
  facet_grid(rows = vars(simtype),
             cols = vars(profile)) +
  scale_y_continuous(labels = scales::percent) +
  theme_custom

# plot diff
p_qreg_diff <- effs_qual_reg %>% 
  filter(rvt == "diff") %>% 
  ggplot(aes(x, predicted, color = nat_haz)) +
  geom_hline(yintercept = 0, color = "darkgrey") +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width    = 0.25,
                position = position_dodge(width = 0.2)) +
  labs(x     = expression(Q[reg]),
       y     = "\u0394 PQ",
       color = "Natural\nhazard") +
  facet_grid(rows = vars(simtype),
             cols = vars(profile)) +
  theme_custom

# combine plots
p_qreg_list <- list()
p_qreg_list[[1]] <- p_qreg_abs + 
  theme(legend.position = "none")
p_qreg_list[[2]] <- p_qreg_diff +
  theme(legend.position = "none")
p_qreg_list[[3]] <- get_legend(p_qreg_abs)

gg_pqreg <- plot_grid(plotlist = p_qreg_list,
                    nrow = 1,
                    rel_widths = c(1, 1, 0.2),
                    labels = c("a", "b", ""))

ggsave(filename = str_c(folder_out, "betareg_full_qreg.jpg"),
       plot     = gg_pqreg,
       width    = 16,
       height   = 8,
       units    = "cm",
       scale    = 1.5)


# plot mgm ----------------------------------------------------------------
plot_effects_comb_mgm_qreg(pr = "MP", rvt = "abs")
plot_effects_comb_mgm_qreg(pr = "IP", rvt = "abs")
plot_effects_comb_mgm_qreg(pr = "MP", rvt = "diff")
plot_effects_comb_mgm_qreg(pr = "IP", rvt = "diff")


