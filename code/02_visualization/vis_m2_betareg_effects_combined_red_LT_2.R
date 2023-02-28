## visualize effects of betareg-models
# red LT combined
# variant 2 (only interval & intensity)
## 28.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(ggeffects)
library(cowplot)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v2/"
folder_out <- "results/vis_models_v2/betareg/vis_betareg_red_effects_combined_2/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in_d, "analysis_data_transf.rds")) %>% 
  filter(q_reg2 == "normal")

# subset data and save as new objects
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
    filter(stratum == pull(vals, stratum))
  
  assign(object_name, object_content)
}

rm(object_name, object_content)


# load models -------------------------------------------------------------
models <- read_rds(str_c(folder_in_m, "models_red_LT.rds"))
model_combinations <- read_rds(str_c(folder_in_m, "model_combinations_red_LT.rds")) %>% 
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

plot_effects_comb <- function(st, rvt, stra,
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
    filter(resp_var_type == rvt) %>% 
    filter(stratum == stra)
  
  # calculate effects
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
  
  effs_all <- bind_rows(effs_mgm_intint1,
                        effs_mgm_intint2)
  
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
  p_title_plot <- plot_profilename(str_c(st, rvt, stra, sep = " "))
  p_title_MP <- plot_profilename("Minimal Profile")
  p_title_IP <- plot_profilename("Ideal Profile")
  p_title_intint1 <- plot_predictorname("Interval * Intensity")
  p_title_intint2 <- plot_predictorname("Intensity * Interval")
  
  
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
  plots <- vector(mode = "list", length = 16)
  
  plots[[2]] <- p_title_plot
  plots[[6]] <- p_title_MP
  plots[[7]] <- p_title_IP
  plots[[9]] <- p_title_intint1
  plots[[10]] <- p_intint1_MP
  plots[[11]] <- p_intint1_IP
  plots[[12]] <- p_legend_nathaz2
  plots[[13]] <- p_title_intint2
  plots[[14]] <- p_intint2_MP
  plots[[15]] <- p_intint2_IP
  
  gg_out <- plot_grid(plotlist    = plots,
                      align       = "h",
                      axis        = "lr",
                      ncol        = 4,
                      rel_heights = c(0.25, 0.25, 2, 2),
                      rel_widths  = c(0.15, 1, 1, 0.4))
  
  ggsave(filename = str_c(f_out,
                          st, "_",
                          rvt, "_",
                          stra, ".jpg"),
         plot     = gg_out,
         width    = 20,
         height   = 22,
         units    = "cm",
         scale    = 1)
}


# plot --------------------------------------------------------------------
plot_combinations <- expand_grid(simtype       = "LT",
                                 resp_var_type = c("abs", "diff"),
                                 stratum       = c("UM", "HM", "SA"))

for (i in 1:nrow(plot_combinations)) {
  cat("\r", i, "/", nrow(plot_combinations))
  vals <- plot_combinations[i, ]
  plot_effects_comb(st    = pull(vals, simtype),
                    rvt   = pull(vals, resp_var_type),
                    stra  = pull(vals, stratum))
}
