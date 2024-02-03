## visualization relaxed assessment raw data
## sha_y_MP_met & sha_y_IP_met
## reduced data
## 3.2.24, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(cowplot)

folder_in <- "data/raw/nais_analysis_data/"
folder_out <- "results/vis_raw_data_relaxed/vis_red_alt/"


# functions ---------------------------------------------------------------
prep_data_red <- function (orig) {
  res_mgm <- orig %>% 
    filter(mgm != "NOM")
  
  res_nom <- orig %>% 
    filter(mgm == "NOM")
  
  res_nom_mgm <- vector(mode = "list", length = nrow(res_nom))
  mgms <- unique(res_mgm$mgm)
  for (i in seq_along(res_nom_mgm)) {
    res_nom_mgm[[i]] <- res_nom[i, ] %>% 
      slice(rep(1, 6)) %>% 
      mutate(mgm = mgms)
  }
  res_nom_mgm <- bind_rows(res_nom_mgm)
  
  res <- bind_rows(res_mgm, res_nom_mgm) %>% 
    filter(q_reg > 1) %>% 
    mutate(stratum = factor(stratum, levels = c("UM", "HM", "SA")),
           mgm     = factor(mgm, levels = c("STS", "GRS1", "GRS2", "CAB1", "CAB2", "SC")),
           q_site2 = case_when(q_site == 3 ~ "medium",
                               TRUE        ~ "good"),
           q_site2 = factor(q_site2, levels = c("good", "medium"))) %>% 
    select(stratum, q_site2, init, mgm, mgm_interval, mgm_intensity,
           nat_haz, sha_y_MP_met, sha_y_IP_met)
  
  return(res)
}


plot_title_horiz <- function(name) {
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

plot_title_vert <- function(name) {
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


plot_raw <- function(sel_stra, sel_init, sel_nathaz, sel_q_site, sel_profile,
                     dat = analysis_red) {
  
  # response variable
  if (sel_profile == "MP") {
    sel_respvar <- "sha_y_MP_met"
  } else if (sel_profile == "IP") {
    sel_respvar <- "sha_y_IP_met"
  }
  
  # reduce data
  dat_red <- dat %>% 
    filter(stratum == sel_stra) %>% 
    filter(init == sel_init) %>% 
    filter(nat_haz == sel_nathaz) %>% 
    filter(q_site2 == sel_q_site)
  
  # calculate midpoint
  mp <- dat_red %>% 
    filter(mgm_interval == 0) %>% 
    pull(sel_respvar) %>% 
    unique()
  
  # extract best option
  dat_red_best <- dat_red %>% 
    filter(.data[[sel_respvar]] == max(.data[[sel_respvar]]))
  
  # plot
  gg_out <- ggplot(dat_red, aes(x     = mgm_interval,
                      y     = mgm_intensity,
                      color = .data[[sel_respvar]])) +
    geom_point(data = dat_red_best, color = "green", size = 4) +
    geom_point(size = 3) +
    scale_color_gradient2(low      = "#ff0000",
                          mid      = "#B2BEB5",
                          high     = "#004a8d",
                          midpoint = mp) +
    facet_wrap(~ mgm) +
    labs(x     = "Interval (y)",
         y     = "Intensity (%)",
         color = NULL) +
    coord_equal() +
    theme_bw() +
    theme(panel.grid = element_blank())
  
  return(gg_out)
}


plot_raw_comb <- function(stratum_sel, init_sel, nathaz_sel,
                          f_out = folder_out) {
  plots <- vector(mode = "list", length = 9)
  plots[[1]] <- plot_title_horiz(name = nathaz_sel)
  plots[[2]] <- plot_title_horiz(name = expression(bold(Q[site]*": good")))
  plots[[3]] <- plot_title_horiz(name = expression(bold(Q[site]*": medium")))
  plots[[4]] <- plot_title_vert(name = "Minimal profile")
  plots[[5]] <- plot_raw(sel_stra    = stratum_sel,
                         sel_init    = init_sel,
                         sel_nathaz  = nathaz_sel,
                         sel_q_site  = "good",
                         sel_profile = "MP")
  
  plots[[6]] <- plot_raw(sel_stra    = stratum_sel,
                         sel_init    = init_sel,
                         sel_nathaz  = nathaz_sel,
                         sel_q_site  = "medium",
                         sel_profile = "MP")
  plots[[7]] <- plot_title_vert(name = "Ideal profile")
  plots[[8]] <- plot_raw(sel_stra    = stratum_sel,
                         sel_init    = init_sel,
                         sel_nathaz  = nathaz_sel,
                         sel_q_site  = "good",
                         sel_profile = "IP")
  
  plots[[9]] <- plot_raw(sel_stra    = stratum_sel,
                         sel_init    = init_sel,
                         sel_nathaz  = nathaz_sel,
                         sel_q_site  = "medium",
                         sel_profile = "IP")
  
  gg_comb <- plot_grid(plotlist = plots,
            align    = "v",
            rel_widths = c(0.1, 1, 1),
            rel_heights = c(0.1, 1, 1))
  
  ggsave(filename = str_c(f_out,
                          stratum_sel,
                          "_init", init_sel,
                          "_", nathaz_sel, "_shayP.jpg"),
         plot     = gg_comb,
         width    = 8,
         height   = 6,
         units    = "cm",
         scale    = 2.5)
  }


# load and prepare data ---------------------------------------------------
analysis_LT_orig <- read_rds(str_c(folder_in, "LT_relaxed.rds"))
LT_red <- prep_data_red(analysis_LT_orig)
analysis_ST_orig <- read_rds(str_c(folder_in, "ST_relaxed.rds"))
ST_red <- prep_data_red(analysis_ST_orig)
analysis_red <- bind_rows(ST_red, LT_red)

rm(analysis_LT_orig, analysis_ST_orig, LT_red, ST_red)


# plot --------------------------------------------------------------------
var_comb <- expand_grid(stratum = c("UM", "HM", "SA"),
                        init    = c("1", "2", "3", "LT"),
                        nathaz  = c("A", "LED"))


for (i in 1:nrow(var_comb)) {
  cat("plot ", i, "/", nrow(var_comb), "\n")
  var_sel <- var_comb[i,]
  
  plot_raw_comb(stratum_sel = pull(var_sel, stratum),
                init_sel    = pull(var_sel, init),
                nathaz_sel  = pull(var_sel, nathaz))
  
}


