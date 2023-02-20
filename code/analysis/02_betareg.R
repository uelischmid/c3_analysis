## analyize data with beta regression
## 20.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(broom)
library(effects)

folder_in <- "data/processed/analysis/"
folder_out <- "results/betareg_effects/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in, "analysis_data_transf.rds"))

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



# prepare models ----------------------------------------------------------
model_combinations <- expand_grid(simtype       = c("LT", "ST"),
                                  nat_haz       = c("A", "LED"),
                                  profile       = c("MP", "IP"),
                                  resp_var_type = c("abs", "diff"),
                                  formula_nr    = c(1:3),
                                  link          = c("logit", "cloglog"))

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
  if (pull(vals, simtype) == "LT") {
    m_formula <- str_c(m_resp_var, " ~ ",
                       "stratum + q_site2 * q_reg2 + mgm_type + mgm_interval + mgm_intensity")
  } else {
    m_formula <- str_c(m_resp_var, " ~ ",
                       "stratum + q_site2 * q_reg2 + init +mgm_type + mgm_interval + mgm_intensity")
  }
  
  if (pull(vals, formula_nr) == 2) {
    m_formula <- str_c(m_formula, " + I(mgm_interval^2) + I(mgm_intensity^2)")
  }
  
  if (pull(vals, formula_nr) == 3) {
    m_formula <- str_c(m_formula, " + I(mgm_interval^2) + I(mgm_intensity^2) + mgm_interval * mgm_intensity")
  }
  
  m_formula <- as.formula(m_formula)
  
  # data
  m_data <- str_c("data_", pull(vals, simtype), "_", pull(vals, nat_haz))
  
  # link
  m_link <- pull(vals, link)
  
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


# compare models ----------------------------------------------------------
glances


# visualize effects -------------------------------------------------------
for (i in seq_along(models)) {
  cat("\r ", i, "/", length(models))
  
  vals <- model_combinations[i, ]
  
  jpeg(filename = str_c(folder_out, "y_auto/br_",
                        pull(vals, simtype), "_",
                        pull(vals, nat_haz), "_",
                        pull(vals, profile), "_",
                        pull(vals, resp_var_type), "_",
                        "f", pull(vals, formula_nr), ".jpg"),
       width    = 1000,
       height   = 700)
  plot(allEffects(mod = models[[i]]))
  dev.off()
  
  jpeg(filename = str_c(folder_out, "y_01/br_",
                        pull(vals, simtype), "_",
                        pull(vals, nat_haz), "_",
                        pull(vals, profile), "_",
                        pull(vals, resp_var_type), "_",
                        "f", pull(vals, formula_nr), ".jpg"),
       width    = 1000,
       height   = 700)
  plot(allEffects(mod = models[[i]]),
       axes = list(y = list(lim = c(0, 1),
                            type = "response")))
  dev.off()
}
