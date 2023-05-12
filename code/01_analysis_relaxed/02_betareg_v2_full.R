## analyze data with beta regression
## full models
## relaxed assessments
## 12.5.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(broom)

folder_in <- "data/processed/nais_analysis_data/"
folder_out <- "data/processed/betareg_models_v2_relaxed/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in, "analysis_data_transf_relaxed.rds"))

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
                                  resp_var_type = c("abs", "diff"))

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
                       "stratum + q_site2 * q_reg2 + mgm_type + mgm_interval * mgm_intensity + I(mgm_interval^2) + I(mgm_intensity^2)")
  } else {
    m_formula <- str_c(m_resp_var, " ~ ",
                       "stratum + q_site2 * q_reg2 + init + mgm_type + mgm_interval * mgm_intensity + I(mgm_interval^2) + I(mgm_intensity^2)")
  }
  m_formula <- as.formula(m_formula)
  
  # data
  m_data <- str_c("data_", pull(vals, simtype), "_", pull(vals, nat_haz))
  
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
ggplot(glances, aes(simtype, pseudo.r.squared)) +
  geom_boxplot() +
  facet_wrap(~resp_var_type)


# save model data ---------------------------------------------------------
write_rds(models, str_c(folder_out, "models_full.rds"))
write_rds(glances, str_c(folder_out, "glances_full.rds"))
write_rds(model_combinations, str_c(folder_out, "model_combinations_full.rds"))
