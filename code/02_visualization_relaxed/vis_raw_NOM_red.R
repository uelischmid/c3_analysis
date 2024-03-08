# visualize protective function of NOM-mgm
# relaxed assessment
# 1 graph per stratum
# 15.1.24

# setup -------------------------------------------------------------------
library(tidyverse)
path_in <- "data/raw/nais_analysis_data/"
path_out <- "results/vis_raw_data_relaxed/vis_NOM_red2/"

# load data ---------------------------------------------------------------
LT_NOM_red <- read_rds(str_c(path_in, "LT_relaxed.rds")) %>% 
  select(stratum:nat_haz, starts_with("sha_i")) %>% 
  filter(q_reg > 1) %>% 
  filter(mgm == "NOM") %>% 
  select(stratum, q_site, nat_haz, sha_i_MP_met, sha_i_IP_met)

ST_NOM_red <- read_rds(str_c(path_in, "ST_relaxed.rds")) %>% 
  select(stratum:nat_haz, starts_with("sha_i")) %>% 
  filter(q_reg > 1) %>% 
  filter(mgm == "NOM") %>% 
  select(stratum, init, q_site, nat_haz, sha_i_MP_met, sha_i_IP_met)

NOM_red <- LT_NOM_red %>% 
  mutate(init = "LT") %>% 
  select(stratum, init, everything()) %>% 
  bind_rows(ST_NOM_red)

# transform data ----------------------------------------------------------
LT_NOM_red_t <- LT_NOM_red %>%
  mutate(stratum = factor(stratum, levels = c("UM", "HM", "SA")),
         q_site2 = case_when(q_site == 3 ~ "medium",
                             TRUE        ~ "good"),
         q_site2 = factor(q_site2, levels = c("good", "medium")),
         nat_haz = factor(nat_haz, levels = c("A", "LED"))) %>% 
  select(stratum, q_site2, nat_haz, sha_i_MP_met, sha_i_IP_met) %>% 
  pivot_longer(cols      = sha_i_MP_met:sha_i_IP_met,
               names_to  = "profile",
               values_to = "pv") %>% 
  mutate(profile = str_sub(profile, 7, 8),
         profile = factor(profile, levels = c("MP", "IP")))
  
ST_NOM_red_t <- ST_NOM_red %>%
  mutate(stratum = factor(stratum, levels = c("UM", "HM", "SA")),
         init2   = case_when(init == "1" ~ "young",
                             init == "2" ~ "structured",
                             TRUE        ~ "mature"),
         init2   = factor(init2, levels = c("young", "structured", "mature")),
         q_site2 = case_when(q_site == 3 ~ "medium",
                             TRUE        ~ "good"),
         q_site2 = factor(q_site2, levels = c("good", "medium")),
         nat_haz = factor(nat_haz, levels = c("A", "LED"))) %>% 
  select(stratum, init2, q_site2, nat_haz, sha_i_MP_met, sha_i_IP_met) %>% 
  pivot_longer(cols      = sha_i_MP_met:sha_i_IP_met,
               names_to  = "profile",
               values_to = "pv") %>% 
  mutate(profile = str_sub(profile, 7, 8),
         profile = factor(profile, levels = c("MP", "IP")))


NOM_red_t <- NOM_red %>%
  mutate(stratum = factor(stratum, levels = c("UM", "HM", "SA")),
         init2   = case_when(init == "1" ~ "ST: young",
                             init == "2" ~ "ST: structured",
                             init == "3" ~ "ST: mature",
                             TRUE        ~ "LT"),
         init2   = factor(init2, levels = c("ST: young", "ST: structured", "ST: mature", "LT")),
         q_site2 = case_when(q_site == 3 ~ "medium",
                             TRUE        ~ "good"),
         q_site2 = factor(q_site2, levels = c("good", "medium")),
         nat_haz = factor(nat_haz, levels = c("A", "LED"))) %>% 
  select(stratum, init2, q_site2, nat_haz, sha_i_MP_met, sha_i_IP_met) %>% 
  pivot_longer(cols      = sha_i_MP_met:sha_i_IP_met,
               names_to  = "profile",
               values_to = "pv") %>% 
  mutate(profile = str_sub(profile, 7, 8),
         profile = factor(profile, levels = c("MP", "IP")))

# plot --------------------------------------------------------------------
plot_nom <- function(stra,
                     df = NOM_red_t, p_out = path_out) {
  gg_out <- df %>% 
    filter(stratum == stra) %>% 
    ggplot(aes(profile, pv, color = nat_haz, shape = q_site2)) +
    geom_point(position = position_dodge2(w = 0.5)) +
    facet_grid(cols = vars(init2)) +
    scale_y_continuous(labels = scales::percent,
                       limits = c(0, 1))+
    labs(x     = NULL,
         y     = "Profile met (%)",
         color = "Natural\nhazard",
         shape = "Site\nquality") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  ggsave(filename = str_c(p_out, "NOM_red_", stra, ".jpg"),
         plot     = gg_out,
         width    = 15,
         height   = 5,
         units    = "cm",
         scale    = 1.3)
}

plot_nom("UM")
plot_nom("HM")
plot_nom("SA")
