## visualize marginal effects of betareg-models
## v2 full
## relaxed assessments
## 12.5.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(margins)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v2_relaxed/"
folder_out <- "results/vis_models_v2_relaxed/betareg/vis_betareg_full_marginal_effects/"


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


# function ----------------------------------------------------------------
plot_mageff <- function(m_all, mc_all,
                        st, rv) {
  mc_used <- mc_all %>% 
    filter(simtype == st) %>% 
    filter(resp_var_type == rv)
  
  margs <- vector(mode = "list", length = nrow(mc_used))
  
  for (i in 1:nrow(mc_used)) {
    vals <- mc_used[i, ]
    
    margs[[i]] <- margins(m_all[[pull(vals, rowid)]], change = "minmax") %>% 
      summary() %>% 
      as_tibble() %>% 
      mutate(nat_haz = pull(vals, nat_haz),
             profile = pull(vals, profile))
  }
  
  margs <- bind_rows(margs) %>% 
  mutate(Predictor = factor,
         Predictor = str_replace_all(Predictor,
                                     c("stratumHM"      = "Stratum: HM",
                                       "stratumSA"      = "Stratum: SA",
                                       "q_site2medium"  = "Qsite: medium",
                                       "q_reg2hindered" = "Qreg: hindered",
                                       "mgm_typeSC"     = "Type: SC",
                                       "mgm_typeGRS2"   = "Type: GRS2",
                                       "mgm_typeGRS1"   = "Type: GRS1",
                                       "mgm_typeCAB2"   = "Type: CAB2",
                                       "mgm_typeCAB1"   = "Type: CAB1",
                                       "mgm_interval"   = "Interval",
                                       "mgm_intensity"  = "Intensity",
                                       "init2"          = "Init: 2",
                                       "init3"          = "Init: 3")),
         Predictor = factor(Predictor, levels = c("Stratum: HM", "Stratum: SA",
                                                  "Init: 2", "Init: 3",
                                                  "Qsite: medium", "Qreg: hindered",
                                                  "Type: GRS1", "Type: GRS2", "Type: CAB1", "Type: CAB2", "Type: SC",
                                                  "Interval", "Intensity")),
         Predictor = fct_rev(Predictor),
         profile   = factor(profile, levels = c("MP", "IP")))
  
  gg_p <- ggplot(margs, aes(Predictor, AME,
                            color    = nat_haz,
                            shape    = profile,
                            linetype = profile)) +
    geom_point(position = position_dodge(width = 0.7)) +
    geom_errorbar(aes(ymin = lower,
                      ymax = upper),
                  position = position_dodge(width = 0.7),
                  width    = 0.25) +
    geom_hline(yintercept = 0) +
    labs(y        = "Average marginal effect",
         color    = "Natural\nhazard",
         shape    = "Profile",
         linetype = "Profile") +
    coord_flip() +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  if (st == "LT") {
    gg_p <- gg_p +
      geom_vline(xintercept = c(1.5, 2.5, 7.5, 8.5, 9.5),
                 linetype   = 2)
  } else if (st == "ST") {
    gg_p <- gg_p +
      geom_vline(xintercept = c(1.5, 2.5, 7.5, 8.5, 9.5, 11.5),
                 linetype   = 2)
  }
  
  return(gg_p)
}


# plot --------------------------------------------------------------------
# LT abs
p_LT_abs <- plot_mageff(m_all  = models,
                        mc_all = model_combinations,
                        st     = "LT",
                        rv     = "abs")

ggsave(filename = str_c(folder_out, "LT_abs.jpg"),
       plot     = p_LT_abs,
       width    = 15,
       height   = 15,
       units    = "cm",
       scale    = 1)


# LT diff
p_LT_diff <- plot_mageff(m_all  = models,
                         mc_all = model_combinations,
                         st     = "LT",
                         rv     = "diff")

ggsave(filename = str_c(folder_out, "LT_diff.jpg"),
       plot     = p_LT_diff,
       width    = 15,
       height   = 15,
       units    = "cm",
       scale    = 1)


# ST abs
p_ST_abs <- plot_mageff(m_all  = models,
                        mc_all = model_combinations,
                        st     = "ST",
                        rv     = "abs")

ggsave(filename = str_c(folder_out, "ST_abs.jpg"),
       plot     = p_ST_abs,
       width    = 15,
       height   = 15,
       units    = "cm",
       scale    = 1)


# ST diff
p_ST_diff <- plot_mageff(m_all  = models,
                         mc_all = model_combinations,
                         st     = "ST",
                         rv     = "diff")

ggsave(filename = str_c(folder_out, "ST_diff.jpg"),
       plot     = p_ST_diff,
       width    = 15,
       height   = 15,
       units    = "cm",
       scale    = 1)
