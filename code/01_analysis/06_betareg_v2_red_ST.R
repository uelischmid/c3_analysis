## analyze data with beta regression
## models with reduced data: ST
## 24.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(broom)

folder_in <- "data/processed/nais_analysis_data/"
folder_out <- "data/processed/betareg_models_v2/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in, "analysis_data_transf.rds")) %>% 
  filter(q_reg2 == "normal")

# subset data and save as new objects
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

# prepare models ----------------------------------------------------------
model_combinations <- expand_grid(simtype       = "ST",
                                  nat_haz       = c("A", "LED"),
                                  profile       = c("MP", "IP"),
                                  resp_var_type = c("abs", "diff"),
                                  stratum       = c("UM", "HM", "SA"),
                                  init          = c("1", "2", "3"))

models <- vector(mode = "list", length = nrow(model_combinations))
glances <- models


# build models ------------------------------------------------------------
for (i in 1:nrow(model_combinations)) {
  cat("\r ", i, "/", nrow(model_combinations))
  vals <- model_combinations[i, ]
  
  # response variable
  if (pull(vals, resp_var_type) == "abs") {
    if (pull(vals, profile) == "MP") {
      m_resp_var <- "sha_i_MP_met_abs"
    } else {
      m_resp_var <- "sha_i_IP_met_abs"
    }
  } else {
    if (pull(vals, profile) == "MP") {
      m_resp_var <- "sha_i_MP_met_diff_t"
    } else {
      m_resp_var <- "sha_i_IP_met_diff_t"
    }
  }
  
  # formula
  m_formula <- str_c(m_resp_var, " ~ ",
                     "q_site2 + mgm_type + mgm_interval * mgm_intensity + I(mgm_interval^2) + I(mgm_intensity^2)")
  m_formula <- as.formula(m_formula)
  
  # data
  m_data <- str_c("data_",
                  pull(vals, simtype), "_",
                  pull(vals, nat_haz), "_",
                  pull(vals, stratum), "_",
                  pull(vals, init))
  
  # link
  m_link <- "logit"
  
  # build model
  m <- eval(
    bquote(betareg(formula = .(m_formula),
                   data    = .(as.name(m_data)),
                   link    = .(m_link)))
  )
  
  # save model
  models[[i]] <- m
  
  # save glance
  glances[[i]] <- bind_cols(vals, glance(m))
}
cat("\n")

glances <- bind_rows(glances)


# plot r squared ----------------------------------------------------------
ggplot(glances, aes(nat_haz, pseudo.r.squared)) +
  geom_boxplot() +
  facet_wrap(~init)


# save model data ---------------------------------------------------------
write_rds(models, str_c(folder_out, "models_red_ST.rds"))
write_rds(glances, str_c(folder_out, "glances_red_ST.rds"))
write_rds(model_combinations, str_c(folder_out, "model_combinations_red_ST.rds"))
