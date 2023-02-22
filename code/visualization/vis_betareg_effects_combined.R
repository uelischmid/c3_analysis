## visualize effects of betareg-models
## 20.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(ggeffects)
library(cowplot)

folder_in_d <- "data/processed/analysis/"
folder_in_m <- "results/betareg_models/"
folder_out <- "results/vis_betareg_effects_combined/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in_d, "analysis_data_transf.rds"))

data_LT_A <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "A")

data_LT_LED <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "LED")

data_ST_A <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "A")

data_ST_LED <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "LED")


# load models -------------------------------------------------------------
models <- read_rds(str_c(folder_in_m, "models.rds"))
model_combinations <- read_rds(str_c(folder_in_m, "model_combinations.rds")) %>% 
  rowid_to_column()



# functions ---------------------------------------------------------------
get_eff <- function(mod, mc, term) {
  effs <- vector(mode = "list", length = nrow(mc))
  
  for (i in 1:nrow(mc)) {
    vals <- mc[i, ]
    effs[[i]] <- ggeffect(model = mod[[pull(vals, rowid)]],
                          terms = term) %>% 
      as_tibble() %>% 
      mutate(nat_haz = pull(vals, nat_haz),
             profile = pull(vals, profile))
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


plot_effects_comb <- function(st, rvt, fn,
                              yl_MP, yl_IP,
                              m_all, mc_all) {
  # select model combinations
  mc_used <- mc_all %>% 
    filter(link == "logit") %>% 
    filter(simtype == st) %>% 
    filter(resp_var_type == rvt) %>% 
    filter(formula_nr == fn)
  
  # prepare plots
  theme_custom <- theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  if (fn %in% c(1, 2)) {
    plots <- vector(mode = "list", length = 24)
  } else if (fn == 3) {
    plots <- vector(mode = "list", length = 32)
  }
  
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
  
  
  # name plots
  plots[[2]] <- plot_profilename("Minimal Profile")
  plots[[3]] <- plot_profilename("Ideal Profile")
  plots[[5]] <- plot_predictorname("Stratum")
  plots[[9]] <- plot_predictorname("Quality")
  plots[[13]] <- plot_predictorname("Intervention\ntype")
  plots[[17]] <- plot_predictorname("Intervention\ninterval")
  plots[[21]] <- plot_predictorname("Intervention\nintensity")
  if (fn == 3) {
    plots[[25]] <- plot_predictorname("Intensity*\nInterval")
    plots[[29]] <- plot_predictorname("Interval*\nIntensity")
  }
  
  # plot stratum
  effs_stratum <- get_eff(mod  = m_all,
                          mc   = mc_used,
                          term = "stratum")
  
  gg_MP_temp <- effs_stratum %>% 
    filter(profile == "MP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x     = "Stratum",
         y     = "Profile met (%)",
         color = "Natural\nhazard") +
    scale_y_continuous(limits = yl_MP,
                       breaks = yb_MP,
                       labels = scales::percent) +
    theme_custom +
    theme(legend.justification = "top")
  
  plots[[6]] <- gg_MP_temp +
    theme(legend.position = "none")
  
  plots[[7]] <- effs_stratum %>% 
    filter(profile == "IP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x = "Stratum",
         y = " ") +
    scale_y_continuous(limits = yl_IP,
                       breaks = yb_IP,
                       labels = scales::percent) +
    theme_custom +
    theme(legend.position = "none")
  
  plots[[8]] <- get_legend(gg_MP_temp)
  
  # plot qualities (x = site)
  effs_quality <- get_eff(mod  = m_all,
                          mc   = mc_used,
                          term = c("q_site2", "q_reg2")) %>% 
    rename(Qreg = group)
  
  gg_MP_temp <- effs_quality %>% 
    filter(profile == "MP") %>% 
    ggplot(aes(x, predicted, color = nat_haz, shape = Qreg, linetype = Qreg)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x        = expression(Q[site]),
         y        = "Profile met (%)",
         color    = "Natural\nhazard",
         shape    = expression(Q[reg]),
         linetype = expression(Q[reg])) +
    scale_y_continuous(limits = yl_MP,
                       breaks = yb_MP,
                       labels = scales::percent) +
    theme_custom +
    scale_color_discrete(guide = "none") +
    theme(legend.justification = "top")
  
  plots[[10]] <- gg_MP_temp +
    theme(legend.position = "none")
  
  plots[[11]] <- effs_quality %>% 
    filter(profile == "IP") %>% 
    ggplot(aes(x, predicted, color = nat_haz, shape = Qreg, linetype = Qreg)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x = expression(Q[site]),
         y = " ") +
    scale_y_continuous(limits = yl_IP,
                       breaks = yb_IP,
                       labels = scales::percent) +
    theme_custom +
    theme(legend.position = "none")
  
  plots[[12]] <- get_legend(gg_MP_temp)
  
  # plot mgm type
  effs_mgm_type <- get_eff(mod = m_all,
                           mc = mc_used,
                           term = "mgm_type")
  
  gg_MP_temp <- effs_mgm_type %>% 
    filter(profile == "MP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x     = "Type",
         y     = "Profile met (%)",
         color = "Natural\nhazard") +
    scale_y_continuous(limits = yl_MP,
                       breaks = yb_MP,
                       labels = scales::percent) +
    theme_custom +
    theme(legend.justification = "top")
  
  plots[[14]] <- gg_MP_temp +
    theme(legend.position = "none")
  
  plots[[15]] <- effs_mgm_type %>% 
    filter(profile == "IP") %>% 
    ggplot(aes(x, predicted, color = nat_haz)) +
    geom_point(position = position_dodge(width = 0.2)) +
    geom_errorbar(aes(ymin = conf.low,
                      ymax = conf.high),
                  width    = 0.25,
                  position = position_dodge(width = 0.2)) +
    labs(x = "Type",
         y = " ") +
    scale_y_continuous(limits = yl_IP,
                       breaks = yb_IP,
                       labels = scales::percent) +
    theme_custom +
    theme(legend.position = "none")
  
  # plots[[16]] <- get_legend(gg_MP_temp)
  
  
  # plot mgm interval
  effs_mgm_interval <- get_eff(mod  = m_all,
                               mc   = mc_used,
                               term = "mgm_interval [10:40]")
  
  gg_MP_temp <- effs_mgm_interval %>% 
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
    labs(x     = "Interval (y)",
         y     = "Profile met (%)",
         color = "Natural\nhazard",
         fill  = "Natural\nhazard") +
    scale_y_continuous(limits = yl_MP,
                       breaks = yb_MP,
                       labels = scales::percent) +
    theme_custom +
    theme(legend.justification = "top")
  
  plots[[18]] <- gg_MP_temp +
    theme(legend.position = "none")
  
  plots[[19]] <- effs_mgm_interval %>% 
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
    labs(x = "Interval (y)",
         y = " ") +
    scale_y_continuous(limits = yl_IP,
                       breaks = yb_IP,
                       labels = scales::percent) +
    theme_custom +
    theme(legend.position = "none")
  
  plots[[20]] <- get_legend(gg_MP_temp)
  
  
  # plot mgm_intensity
  effs_mgm_intensity <- get_eff(mod  = m_all,
                                mc   = mc_used,
                                term = "mgm_intensity [10:40]")
  
  gg_MP_temp <- effs_mgm_intensity %>% 
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
    labs(x     = "Intensity (%)",
         y     = "Profile met (%)",
         color = "Natural\nhazard",
         fill  = "Natural\nhazard") +
    scale_y_continuous(limits = yl_MP,
                       breaks = yb_MP,
                       labels = scales::percent) +
    theme_custom +
    theme(legend.justification = "top")
  
  plots[[22]] <- gg_MP_temp +
    theme(legend.position = "none")
  
  plots[[23]] <- effs_mgm_intensity %>% 
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
    labs(x = "Intensity (%)",
         y = " ") +
    scale_y_continuous(limits = yl_IP,
                       breaks = yb_IP,
                       labels = scales::percent) +
    theme_custom +
    theme(legend.position = "none")
  
  # plots[[24]] <- get_legend(gg_MP_temp)
  
  if (fn == 3) {
    # plot int-int-interaction (x = intensity)
    effs_mgm_interaction1 <- get_eff(mod  = m_all,
                                     mc   = mc_used,
                                     term = c("mgm_intensity [10:40]",
                                              "mgm_interval [10, 20, 30, 40]")) %>% 
      rename(Interval = group)
    
    gg_MP_temp <- effs_mgm_interaction1 %>% 
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
           y     = "Profile met (%)",
           color = "Natural\nhazard",
           fill  = "Natural\nhazard") +
      scale_y_continuous(limits = yl_MP,
                         breaks = yb_MP,
                         labels = scales::percent) +
      theme_custom +
      theme(legend.justification = "top")
    
    plots[[26]] <- gg_MP_temp +
      theme(legend.position = "none")
    
    plots[[27]] <- effs_mgm_interaction1 %>% 
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
      labs(x     = "Intensity (%)",
           y     = "Profile met (%)",
           color = "Natural\nhazard",
           fill  = "Natural\nhazard") +
      scale_y_continuous(limits = yl_IP,
                         breaks = yb_IP,
                         labels = scales::percent) +
      theme_custom +
      theme(legend.position = "none")
      
    # plots[[28]] <- get_legend(gg_MP_temp)
    
    
    # plot int-int-interaction (x = interval)
    effs_mgm_interaction2 <- get_eff(mod  = m_all,
                                     mc   = mc_used,
                                     term = c("mgm_interval [10:40]",
                                              "mgm_intensity [10, 20, 30, 40]")) %>% 
      rename(Intensity = group)
    
    gg_MP_temp <- effs_mgm_interaction2 %>% 
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
           y     = "Profile met (%)",
           color = "Natural\nhazard",
           fill  = "Natural\nhazard") +
      scale_y_continuous(limits = yl_MP,
                         breaks = yb_MP,
                         labels = scales::percent) +
      theme_custom +
      theme(legend.justification = "top")
    
    plots[[30]] <- gg_MP_temp +
      theme(legend.position = "none")
    
    plots[[31]] <- effs_mgm_interaction2 %>% 
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
      labs(x     = "Interval (y)",
           y     = "Profile met (%)",
           color = "Natural\nhazard",
           fill  = "Natural\nhazard") +
      scale_y_continuous(limits = yl_IP,
                         breaks = yb_IP,
                         labels = scales::percent) +
      theme_custom +
      theme(legend.position = "none")
    
    # plots[[32]] <- get_legend(gg_MP_temp)
  }
  
  
  # plot together
  if (fn %in% c(1, 2)) {
    gg_out <- plot_grid(plotlist    = plots,
                        align       = "h",
                        axis        = "lr",
                        ncol        = 4,
                        rel_heights = c(0.25, rep(1, 5)),
                        rel_widths  = c(0.15, 1, 1, 0.4))
  } else if (fn == 3) {
    gg_out <- plot_grid(plotlist    = plots,
                        align       = "h",
                        axis        = "lr",
                        ncol        = 4,
                        rel_heights = c(0.25, rep(1, 7)),
                        rel_widths  = c(0.15, 1, 1, 0.4))
  }
  
  return(gg_out)
}


# plot --------------------------------------------------------------------
# LT abs
p <- plot_effects_comb(st     = "LT",
                       rvt    = "abs",
                       fn     = 1,
                       yl_MP  = c(0.3, 0.6),
                       yl_IP  = c(0.2, 0.4),
                       m_all  = models,
                       mc_all = model_combinations)

ggsave(filename = str_c(folder_out, "LT_abs_f1.jpg"),
       plot     = p,
       width    = 20,
       height   = 30,
       units    = "cm",
       scale    = 1)


# LT diff
p <- plot_effects_comb(st     = "LT",
                       rvt    = "diff",
                       fn     = 1,
                       yl_MP  = c(0.4, 0.55),
                       yl_IP  = c(0.45, 0.6),
                       m_all  = models,
                       mc_all = model_combinations)

ggsave(filename = str_c(folder_out, "LT_diff_f1.jpg"),
       plot     = p,
       width    = 20,
       height   = 30,
       units    = "cm",
       scale    = 1)


# ST abs
p <- plot_effects_comb(st     = "ST",
                       rvt    = "abs",
                       fn     = 2,
                       yl_MP  = c(0.3, 0.65),
                       yl_IP  = c(0.15, 0.45),
                       m_all  = models,
                       mc_all = model_combinations)

ggsave(filename = str_c(folder_out, "ST_abs_f2.jpg"),
       plot     = p,
       width    = 20,
       height   = 30,
       units    = "cm",
       scale    = 1)


# ST diff
p <- plot_effects_comb(st     = "ST",
                       rvt    = "diff",
                       fn     = 2,
                       yl_MP  = c(0.42, 0.49),
                       yl_IP  = c(0.44, 0.5),
                       m_all  = models,
                       mc_all = model_combinations)

ggsave(filename = str_c(folder_out, "ST_diff_f2.jpg"),
       plot     = p,
       width    = 20,
       height   = 30,
       units    = "cm",
       scale    = 1)

