## visualize marginal effects of betareg-models
## 22.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(margins)

folder_in_d <- "data/processed/analysis/"
folder_in_m <- "data/processed/betareg_models_v1/"
folder_out <- "results/vis_models_v1/betareg/vis_betareg_marginal_effects_combined/"


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
  filter(nat_haz == "A") %>% 
  mutate(init = fct_drop(init))

data_ST_LED <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "LED") %>% 
  mutate(init = fct_drop(init))


# load models -------------------------------------------------------------
models <- read_rds(str_c(folder_in_m, "models.rds"))
model_combinations <- read_rds(str_c(folder_in_m, "model_combinations.rds")) %>% 
  rowid_to_column()


# function ----------------------------------------------------------------
plot_mageff <- function(m_all, mc_all,
                        st, rv, fn) {
  mc_used <- mc_all %>% 
    filter(link == "logit") %>% 
    filter(simtype == st) %>% 
    filter(resp_var_type == rv) %>% 
    filter(formula_nr == fn)
  
  margs <- vector(mode = "list", length = nrow(mc_used))
  
  for (i in 1:nrow(mc_used)) {
    vals <- mc_used[i, ]
    
    margs[[i]] <- margins(m_all[[pull(vals, rowid)]], change = "minmax") %>% 
      summary() %>% 
      as_tibble() %>% 
      mutate(nat_haz = pull(vals, nat_haz),
             profile = pull(vals, profile))
  }
  
  margs <- bind_rows(margs)%>% 
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
  
  gg_p <- ggplot(margs, aes(Predictor, AME, color = nat_haz, shape = profile, linetype = profile)) +
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
  
  if (fn == 1) {
    gg_p <- gg_p +
      geom_vline(xintercept = c(1.5, 2.5, 7.5, 8.5, 9.5),
                 linetype   = 2)
  } else if (fn == 2) {
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
                        rv     = "abs",
                        fn     = 1)

ggsave(filename = str_c(folder_out, "LT_abs_f1.jpg"),
       plot     = p_LT_abs,
       width    = 15,
       height   = 15,
       units    = "cm",
       scale    = 1)


# LT diff
p_LT_diff <- plot_mageff(m_all  = models,
                         mc_all = model_combinations,
                         st     = "LT",
                         rv     = "diff",
                         fn     = 1)

ggsave(filename = str_c(folder_out, "LT_diff_f1.jpg"),
       plot     = p_LT_diff,
       width    = 15,
       height   = 15,
       units    = "cm",
       scale    = 1)


# ST abs
p_ST_abs <- plot_mageff(m_all  = models,
                        mc_all = model_combinations,
                        st     = "ST",
                        rv     = "abs",
                        fn     = 2)

ggsave(filename = str_c(folder_out, "ST_abs_f2.jpg"),
       plot     = p_ST_abs,
       width    = 15,
       height   = 15,
       units    = "cm",
       scale    = 1)


# ST diff
p_ST_diff <- plot_mageff(m_all  = models,
                         mc_all = model_combinations,
                         st     = "ST",
                         rv     = "diff",
                         fn     = 2)

ggsave(filename = str_c(folder_out, "ST_diff_f2.jpg"),
       plot     = p_ST_diff,
       width    = 15,
       height   = 15,
       units    = "cm",
       scale    = 1)
