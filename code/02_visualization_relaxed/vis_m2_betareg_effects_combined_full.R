## visualize effects of betareg-models
## full models
## relaxed assessments
## 2.11.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(ggeffects)
library(cowplot)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v2_relaxed/"
folder_out <- "results/vis_models_v2_relaxed/betareg/vis_betareg_full_effects_combined/"


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
             profile = pull(vals, profile))
    
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


plot_effects_comb_fix <- function(st, rvt,
                                  m_all  = models,
                                  mc_all = model_combinations,
                                  f_out  = folder_out) {
  # function for determining scale limits
  get_mm <- function(df) {
    mm <- list()
    
    MP_min <- df %>% 
      filter(profile == "MP") %>% 
      pull(conf.low) %>% 
      min()
    mm$MP$min <- floor(MP_min * 20) / 20
    
    MP_max <- df %>% 
      filter(profile == "MP") %>% 
      pull(conf.high) %>% 
      max()
    mm$MP$max <- ceiling(MP_max * 20) / 20
    
    IP_min <- df %>% 
      filter(profile == "IP") %>% 
      pull(conf.low) %>% 
      min()
    mm$IP$min <- floor(IP_min * 20) /20
    
    IP_max <- df %>% 
      filter(profile == "IP") %>% 
      pull(conf.high) %>% 
      max()
    mm$IP$max <- ceiling(IP_max * 20) / 20
    
    return(mm)
  }
  
  # select model combinations
  mc_used <- mc_all %>% 
    filter(simtype == st) %>% 
    filter(resp_var_type == rvt)
  
  # calculate effects
  effs_stratum <- get_eff(mod  = m_all,
                          mc   = mc_used,
                          term = "stratum")
  effs_qual_site <- get_eff(mod  = m_all,
                            mc   = mc_used,
                            term = "q_site2")
  effs_qual_reg <- get_eff(mod  = m_all,
                           mc   = mc_used,
                           term = "q_reg2")
  if (st == "ST") {
    effs_init <- get_eff(mod  = m_all,
                         mc   = mc_used,
                         term = "init")
  }
  
  
  # prepare plots
  theme_custom <- theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  if (st == "ST") {
    effs_all <- bind_rows(effs_stratum,
                        effs_qual_site,
                        effs_qual_reg,
                        effs_init)
  } else {
    effs_all <- bind_rows(effs_stratum,
                          effs_qual_site,
                          effs_qual_reg)
  }
  mm_all <- get_mm(effs_all)
  
  yl_MP <- c(mm_all$MP$min, mm_all$MP$max)
  yl_IP <- c(mm_all$IP$min, mm_all$IP$max)
  
  if (yl_MP[2] - yl_MP[1] < 0.1) {
    yb_MP <- seq(yl_MP[1], yl_MP[2], 0.02)
  } else if (yl_MP[2] - yl_MP[1] < 0.3) {
    yb_MP <- seq(yl_MP[1], yl_MP[2], 0.05)
  } else {
    yb_MP <- seq(yl_MP[1], yl_MP[2], 0.1)
  }
  
  if (yl_IP[2] - yl_IP[1] < 0.1) {
    yb_IP <- seq(yl_IP[1], yl_IP[2], 0.02)
  } else if (yl_IP[2] - yl_IP[1] < 0.3) {
    yb_IP <- seq(yl_IP[1], yl_IP[2], 0.05)
  } else {
    yb_IP <- seq(yl_IP[1], yl_IP[2], 0.1)
  }
  
  if (rvt == "abs") {
    label_y <- "Profile met (%)"
    scale_y_custom_MP <- scale_y_continuous(limits = yl_MP,
                                            breaks = yb_MP,
                                            labels = scales::percent)
    scale_y_custom_IP <- scale_y_continuous(limits = yl_IP,
                                            breaks = yb_IP,
                                            labels = scales::percent)
  } else {
    label_y <- "PQ difference"
    yb_MP <- round(yb_MP, 1)
    scale_y_custom_MP <- scale_y_continuous(limits = yl_MP,
                                            breaks = yb_MP)
    yb_IP <- round(yb_IP, 1)
    scale_y_custom_IP <- scale_y_continuous(limits = yl_IP,
                                            breaks = yb_IP)
  }
  
  # name plots
  p_title_plot <- plot_profilename(str_c(st, rvt, "(site variables)", sep = " "))
  p_title_MP <- plot_profilename("Minimal Profile")
  p_title_IP <- plot_profilename("Ideal Profile")
  p_title_stratum <- plot_predictorname("Elevational zone")
  p_title_qual_site <- plot_predictorname("Site quality")
  p_title_qual_reg <- plot_predictorname("Regeneration quality")
  p_title_init <- plot_predictorname("Initialization Stand")
  
  
  # plot stratum
  gg_MP_temp <- effs_stratum %>%
    filter(profile == "MP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x        = "Elevational zone",
         y        = label_y,
         color    = "Natural\nhazard") +
    scale_y_custom_MP +
    theme_custom +
    theme(legend.justification = "top")
  
  p_stratum_MP <- gg_MP_temp +
    theme(legend.position = "none")
  
  p_stratum_IP <- effs_stratum %>% 
    filter(profile == "IP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x = "Elevational zone",
         y = " ") +
    scale_y_custom_IP +
    theme_custom +
    theme(legend.position = "none")
  
  p_legend_nathaz1 <- get_legend(gg_MP_temp)
  
  
  # plot site quality
  gg_MP_temp <- effs_qual_site %>%
    filter(profile == "MP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x        = expression(Q[site]),
         y        = label_y,
         color    = "Natural\nhazard") +
    scale_y_custom_MP +
    theme_custom +
    theme(legend.justification = "top")
  
  p_qual_site_MP <- gg_MP_temp +
    theme(legend.position = "none")
  
  p_qual_site_IP <- effs_qual_site %>% 
    filter(profile == "IP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x = expression(Q[site]),
         y = " ") +
    scale_y_custom_IP +
    theme_custom +
    theme(legend.position = "none")
  
  
  # plot regeneration quality
  gg_MP_temp <- effs_qual_reg %>%
    filter(profile == "MP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x        = expression(Q[reg]),
         y        = label_y,
         color    = "Natural\nhazard") +
    scale_y_custom_MP +
    theme_custom +
    theme(legend.justification = "top")
  
  p_qual_reg_MP <- gg_MP_temp +
    theme(legend.position = "none")
  
  p_qual_reg_IP <- effs_qual_reg %>% 
    filter(profile == "IP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x = expression(Q[reg]),
         y = " ") +
    scale_y_custom_IP +
    theme_custom +
    theme(legend.position = "none")
  
  # plot regeneration quality
  if (st == "ST") {
    gg_MP_temp <- effs_init %>%
      filter(profile == "MP") %>% 
      ggplot(aes(x, predicted, color = nat_haz)) +
      geom_point(position = position_dodge(width = 0.2)) +
      geom_errorbar(aes(ymin = conf.low,
                        ymax = conf.high),
                    width    = 0.25,
                    position = position_dodge(width = 0.2)) +
      labs(x        = "Stand type",
           y        = label_y,
           color    = "Natural\nhazard") +
      scale_y_custom_MP +
      theme_custom +
      theme(legend.justification = "top")
    
    p_init_MP <- gg_MP_temp +
      theme(legend.position = "none")
    
    p_init_IP <- effs_init %>% 
      filter(profile == "IP") %>% 
      ggplot(aes(x, predicted, color = nat_haz)) +
      geom_point(position = position_dodge(width = 0.2)) +
      geom_errorbar(aes(ymin = conf.low,
                        ymax = conf.high),
                    width    = 0.25,
                    position = position_dodge(width = 0.2)) +
      labs(x = "Stand type",
           y = " ") +
      scale_y_custom_IP +
      theme_custom +
      theme(legend.position = "none")
    
  } else {
    p_init_MP <- ggplot() +
      theme_void()
    p_init_IP <- p_init_MP
  }
  
   
  # plot together
  plots <- vector(mode = "list", length = 24)
  
  plots[[2]] <- p_title_plot
  plots[[6]] <- p_title_MP
  plots[[7]] <- p_title_IP
  plots[[9]] <- p_title_stratum
  plots[[10]] <- p_stratum_MP
  plots[[11]] <- p_stratum_IP
  plots[[12]] <- p_legend_nathaz1
  plots[[13]] <- p_title_qual_site
  plots[[14]] <- p_qual_site_MP
  plots[[15]] <- p_qual_site_IP
  plots[[17]] <- p_title_qual_reg
  plots[[18]] <- p_qual_reg_MP
  plots[[19]] <- p_qual_reg_IP
  plots[[21]] <- p_title_init
  plots[[22]] <- p_init_MP
  plots[[23]] <- p_init_IP
  
  gg_out <- plot_grid(plotlist    = plots,
                      align       = "h",
                      axis        = "lr",
                      ncol        = 4,
                      rel_heights = c(0.25, 0.25, 1, 1, 1, 1),
                      rel_widths  = c(0.15, 1, 1, 0.4))
  
  ggsave(filename = str_c(f_out,
                          st, "_",
                          rvt, "_sitevar.jpg"),
         plot     = gg_out,
         width    = 20,
         height   = 32,
         units    = "cm",
         scale    = 1)
}


plot_effects_comb_fix(st = "ST",
                      rvt = "diff")


plot_effects_comb_mgm <- function(st, rvt,
                                  m_all  = models,
                                  mc_all = model_combinations,
                                  f_out  = folder_out) {
  # function for determining scale limits
  get_mm <- function(df) {
    mm <- list()
    
    MP_min <- df %>% 
      filter(profile == "MP") %>% 
      pull(conf.low) %>% 
      min()
    mm$MP$min <- floor(MP_min * 20) / 20
    
    MP_max <- df %>% 
      filter(profile == "MP") %>% 
      pull(conf.high) %>% 
      max()
    mm$MP$max <- ceiling(MP_max * 20) / 20
    
    IP_min <- df %>% 
      filter(profile == "IP") %>% 
      pull(conf.low) %>% 
      min()
    mm$IP$min <- floor(IP_min * 20) /20
    
    IP_max <- df %>% 
      filter(profile == "IP") %>% 
      pull(conf.high) %>% 
      max()
    mm$IP$max <- ceiling(IP_max * 20) / 20
    
    return(mm)
  }
  
  # select model combinations
  mc_used <- mc_all %>% 
    filter(simtype == st) %>% 
    filter(resp_var_type == rvt)
  
  # calculate effects
  effs_mgm_type <- get_eff(mod  = m_all,
                           mc   = mc_used,
                           term = "mgm_type")
  effs_mgm_intint1 <- get_eff(mod  = m_all,
                              mc   = mc_used,
                              term = c("mgm_interval [10:40]",
                                       "mgm_intensity [10, 20, 30, 40]"))
  effs_mgm_intint2 <- get_eff(mod  = m_all,
                              mc   = mc_used,
                              term = c("mgm_intensity [10:40]",
                                       "mgm_interval [10, 20, 30, 40]"))
  
  # prepare plots
  theme_custom <- theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  effs_all <- bind_rows(effs_mgm_type,
                        effs_mgm_intint1 %>% 
                          mutate(x = factor(x)),
                        effs_mgm_intint2 %>% 
                          mutate(x = factor(x)))
  
  mm_all <- get_mm(effs_all)
  
  yl_MP <- c(mm_all$MP$min, mm_all$MP$max)
  yl_IP <- c(mm_all$IP$min, mm_all$IP$max)
  
  if (yl_MP[2] - yl_MP[1] < 0.1) {
    yb_MP <- seq(yl_MP[1], yl_MP[2], 0.02)
  } else if (yl_MP[2] - yl_MP[1] < 0.3) {
    yb_MP <- seq(yl_MP[1], yl_MP[2], 0.05)
  } else {
    yb_MP <- seq(yl_MP[1], yl_MP[2], 0.1)
  }
  
  if (yl_IP[2] - yl_IP[1] < 0.1) {
    yb_IP <- seq(yl_IP[1], yl_IP[2], 0.02)
  } else if (yl_IP[2] - yl_IP[1] < 0.3) {
    yb_IP <- seq(yl_IP[1], yl_IP[2], 0.05)
  } else {
    yb_IP <- seq(yl_IP[1], yl_IP[2], 0.1)
  }
  
  if (rvt == "abs") {
    label_y <- "Profile met (%)"
    scale_y_custom_MP <- scale_y_continuous(limits = yl_MP,
                                            breaks = yb_MP,
                                            labels = scales::percent)
    scale_y_custom_IP <- scale_y_continuous(limits = yl_IP,
                                            breaks = yb_IP,
                                            labels = scales::percent)
  } else {
    label_y <- "PQ difference"
    yb_MP <- round(yb_MP, 1)
    scale_y_custom_MP <- scale_y_continuous(limits = yl_MP,
                                            breaks = yb_MP)
    yb_IP <- round(yb_IP, 1)
    scale_y_custom_IP <- scale_y_continuous(limits = yl_IP,
                                            breaks = yb_IP)
  }
  
  # name plots
  p_title_plot <- plot_profilename(str_c(st, rvt, "(management variables)", sep = " "))
  p_title_MP <- plot_profilename("Minimal Profile")
  p_title_IP <- plot_profilename("Ideal Profile")
  p_title_type <- plot_predictorname("Type")
  p_title_intint1 <- plot_predictorname("Interval * Intensity")
  p_title_intint2 <- plot_predictorname("Intensity * Interval")
  
  
  # plot mgm type
  gg_MP_temp <- effs_mgm_type %>% 
    filter(profile == "MP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x     = "Type",
         y     = label_y,
         color = "Natural\nhazard") +
    scale_y_custom_MP +
    theme_custom +
    theme(legend.justification = "top")
  
  p_type_MP <- gg_MP_temp +
    theme(legend.position = "none")
  
  p_type_IP <- effs_mgm_type %>% 
    filter(profile == "IP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x = "Type",
         y = " ") +
    scale_y_custom_IP +
    theme_custom +
    theme(legend.position = "none")
  
  p_legend_nathaz1 <- get_legend(gg_MP_temp)
  
  
  # plot intensity * interval interaction (x = interval)
  effs_mgm_intint1 <- effs_mgm_intint1 %>% 
    rename(Intensity = group)
  
  
  gg_MP_temp <- effs_mgm_intint1 %>% 
    filter(profile == "MP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_ribbon(aes(x, predicted,
                    ymin = conf.low,
                    ymax = conf.high,
                    fill = nat_haz),
                alpha       = 0.2,
                color       = NA,
                inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(~Intensity,
               labeller = "label_both") +
    labs(x     = "Interval (y)",
         y     = label_y,
         color = "Natural\nhazard",
         fill  = "Natural\nhazard") +
    scale_y_custom_MP +
    theme_custom +
    theme(legend.justification = "top")
  
  p_intint1_MP <- gg_MP_temp +
    theme(legend.position = "none")
  
  p_intint1_IP <- effs_mgm_intint1 %>%
    filter(profile == "IP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_ribbon(aes(x, predicted,
                    ymin = conf.low,
                    ymax = conf.high,
                    fill = nat_haz),
                alpha       = 0.2,
                color       = NA,
                inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(~Intensity,
               labeller = "label_both") +
    labs(x = "Interval (y)",
         y = " ") +
    scale_y_custom_IP +
    theme_custom +
    theme(legend.position = "none")
  
  p_legend_nathaz2 <- get_legend(gg_MP_temp)
  
  
  # plot intensity * interval interaction (x = intensity)
  effs_mgm_intint2 <- effs_mgm_intint2 %>% 
    rename(Interval = group)
  
  gg_MP_temp <- effs_mgm_intint2 %>% 
    filter(profile == "MP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_ribbon(aes(x, predicted,
                    ymin = conf.low,
                    ymax = conf.high,
                    fill = nat_haz),
                alpha       = 0.2,
                color       = NA,
                inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(~Interval,
               labeller = "label_both") +
    labs(x     = "Intensity (%)",
         y     = label_y,
         color = "Natural\nhazard",
         fill  = "Natural\nhazard") +
    scale_y_custom_MP +
    theme_custom +
    theme(legend.justification = "top")
  
  p_intint2_MP <- gg_MP_temp +
    theme(legend.position = "none")
  
  p_intint2_IP <- effs_mgm_intint2 %>%
    filter(profile == "IP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_ribbon(aes(x, predicted,
                    ymin = conf.low,
                    ymax = conf.high,
                    fill = nat_haz),
                alpha       = 0.2,
                color       = NA,
                inherit.aes = FALSE) +
    geom_line() +
    facet_wrap(~Interval,
               labeller = "label_both") +
    labs(x = "Intensity (%)",
         y = " ") +
    scale_y_custom_IP +
    theme_custom +
    theme(legend.position = "none")
  
  
  # plot together
  plots <- vector(mode = "list", length = 20)
  
  plots[[2]] <- p_title_plot
  plots[[6]] <- p_title_MP
  plots[[7]] <- p_title_IP
  plots[[9]] <- p_title_type
  plots[[10]] <- p_type_MP
  plots[[11]] <- p_type_IP
  plots[[12]] <- p_legend_nathaz1
  plots[[13]] <- p_title_intint1
  plots[[14]] <- p_intint1_MP
  plots[[15]] <- p_intint1_IP
  plots[[16]] <- p_legend_nathaz2
  plots[[17]] <- p_title_intint2
  plots[[18]] <- p_intint2_MP
  plots[[19]] <- p_intint2_IP
  
  gg_out <- plot_grid(plotlist    = plots,
                      align       = "h",
                      axis        = "lr",
                      ncol        = 4,
                      rel_heights = c(0.25, 0.25, 1, 2, 2),
                      rel_widths  = c(0.15, 1, 1, 0.4))
  
  ggsave(filename = str_c(f_out,
                          st, "_",
                          rvt, "_mgmvar.jpg"),
         plot     = gg_out,
         width    = 20,
         height   = 32,
         units    = "cm",
         scale    = 1)
}

plot_effects_comb_mgm(st = "LT",
                      rvt = "abs")



# plot --------------------------------------------------------------------
plot_combinations <- expand_grid(simtype       = c("ST", "LT"),
                                 resp_var_type = c("abs", "diff"))

for (i in 1:nrow(plot_combinations)) {
  cat("\r", i, "/", nrow(plot_combinations))
  vals <- plot_combinations[i, ]
  plot_effects_comb_fix(st    = pull(vals, simtype),
                        rvt   = pull(vals, resp_var_type))
  plot_effects_comb_mgm(st    = pull(vals, simtype),
                        rvt   = pull(vals, resp_var_type))
}

