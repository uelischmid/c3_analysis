## visualize marginal effects of betareg-models
## v2 red ST
## variant 2 (only qsite & type)
## 28.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(margins)

folder_in_d <- "data/processed/nais_analysis_data/"
folder_in_m <- "data/processed/betareg_models_v2/"
folder_out <- "results/vis_models_v2/betareg/vis_betareg_red_marginal_effects_2/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in_d, "analysis_data_transf.rds")) %>% 
  filter(q_reg2 == "normal")

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
    filter(init    == pull(vals, init))
  
  assign(object_name, object_content)
}

rm(object_name, object_content)


# load models -------------------------------------------------------------
models <- read_rds(str_c(folder_in_m, "models_red_ST.rds"))
model_combinations <- read_rds(str_c(folder_in_m, "model_combinations_red_ST.rds")) %>% 
  rowid_to_column()


# function ----------------------------------------------------------------
plot_mageff <- function(m_all, mc_all, rv) {
  mc_used <- mc_all %>% 
    filter(resp_var_type == rv)
  
  margs <- vector(mode = "list", length = nrow(mc_used))
  
  for (i in 1:nrow(mc_used)) {
    vals <- mc_used[i, ]
    
    margs[[i]] <- margins(m_all[[pull(vals, rowid)]], change = "minmax") %>% 
      summary() %>% 
      as_tibble() %>% 
      mutate(nat_haz = pull(vals, nat_haz),
             profile = pull(vals, profile),
             stratum = pull(vals, stratum),
             init    = pull(vals, init))
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
           profile   = factor(profile, levels = c("MP", "IP")),
           stratum   = factor(stratum, levels = c("UM", "HM", "SA"))) %>% 
    rename(Stratum = stratum,
           Init    = init) %>% 
    filter(!Predictor %in% c("Interval", "Intensity"))
  
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
    geom_vline(xintercept = c(5.5),
               linetype   = 2) +
    labs(title    = str_c("ST ", rv),
         y        = "Average marginal effect",
         color    = "Natural\nhazard",
         shape    = "Profile",
         linetype = "Profile") +
    coord_flip() +
    facet_grid(rows     = vars(Init),
               cols     = vars(Stratum),
               labeller = "label_both") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  
  return(gg_p)
}


# plot --------------------------------------------------------------------
# abs
p_abs <- plot_mageff(m_all  = models,
                     mc_all = model_combinations,
                     rv     = "abs")

ggsave(filename = str_c(folder_out, "ST_abs.jpg"),
       plot     = p_abs,
       width    = 25,
       height   = 25,
       units    = "cm",
       scale    = 1)


# diff
p_diff <- plot_mageff(m_all  = models,
                      mc_all = model_combinations,
                      rv     = "diff")

ggsave(filename = str_c(folder_out, "ST_diff.jpg"),
       plot     = p_diff,
       width    = 25,
       height   = 25,
       units    = "cm",
       scale    = 1)
