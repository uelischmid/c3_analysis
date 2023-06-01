## visualization of LT
## relaxed assessments
## reduced data, sha_i_MP_met & sha_i_IP_met
## 1.6.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(cowplot)

folder_in <- "data/raw/nais_analysis_data/"
folder_out <- "results/vis_raw_data_relaxed/vis_LT_red/"

source("code/02_visualization_relaxed/vis_functions.R")


# load and prepare data ---------------------------------------------------
analysis_LT_orig <- read_rds(str_c(folder_in, "LT_relaxed.rds"))
LT <- prep_data(analysis_LT_orig)
LT_red <- LT %>% 
  filter(q_reg > 1) %>% 
  select(stratum:nat_haz, sha_i_MP_met, sha_i_IP_met) %>%
  mutate(quality = fct_drop(quality))


# plot --------------------------------------------------------------------
plot_combinations <- expand_grid(stratum = c("UM", "HM", "SA"),
                                 init    = "LT",
                                 nat_haz = c("A", "LED"))

for (i in 1:nrow(plot_combinations)) {
  cat(i, "/", nrow(plot_combinations), "\n")
  
  # plot with different color scales
  gg_4cs <- plot_mgm_4cs_all_red(stratum_sel  = plot_combinations[[i, "stratum"]],
                                 init_sel     = plot_combinations[[i, "init"]],
                                 nathaz_sel   = plot_combinations[[i, "nat_haz"]],
                                 res          = LT_red)
  ggsave(filename = str_c(folder_out, "LT_",
                          plot_combinations[[i, "stratum"]], "_",
                          plot_combinations[[i, "nat_haz"]], "_4cs.jpg"),
         plot     = gg_4cs,
         width    = 15,
         height   = 15,
         units    = "cm",
         scale    = 1.2)
 }

