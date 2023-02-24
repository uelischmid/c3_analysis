## analyze data with beta regression
## 20.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
library(broom)
library(effects)

folder_in <- "data/processed/analysis/"
folder_out <- "data/processed/betareg_models/"


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


# save model data ---------------------------------------------------------
write_rds(models, str_c(folder_out, "models.rds"))
write_rds(glances, str_c(folder_out, "glances.rds"))
write_rds(model_combinations, str_c(folder_out, "model_combinations.rds"))



# load model data ---------------------------------------------------------
models <- read_rds(str_c(folder_out, "models.rds"))
glances <- read_rds(str_c(folder_out, "glances.rds"))


# compare link functions --------------------------------------------------
## pseudo r squared
link_diff_r2 <- glances %>% 
  select(-c(df.null:nobs)) %>% 
  pivot_wider(names_from  = link,
              values_from = pseudo.r.squared) %>% 
  mutate(logit_cloglog = logit - cloglog,
         better_r2     = case_when(logit_cloglog > 0 ~ "logit",
                                   logit_cloglog < 0 ~ "cloglog",
                                   TRUE              ~ "tied"))
link_diff_r2 %>% 
  count(better_r2) # equal numbers

ggplot(link_diff_r2, aes(better_r2, abs(logit_cloglog))) +
  geom_boxplot() # if logit better, larger difference

link_diff_r2 %>% 
  count(simtype, better_r2) # no differentiation

link_diff_r2 %>% 
  count(nat_haz, better_r2) # A: cloglog better; LED: logit better

link_diff_r2 %>% 
  count(profile, better_r2) # no differentiation

link_diff_r2 %>% 
  count(resp_var_type, better_r2) # abs: logit better; diff: cloglog better (larger difference than nat_haz)

link_diff_r2 %>% 
  count(formula_nr, better_r2) # no differentiation

# probably use logit for abs and cloglog for diff

## AIC
link_diff_AIC <- glances %>% 
  select(simtype:link, AIC) %>% 
  pivot_wider(names_from  = link,
              values_from = AIC) %>% 
  mutate(logit_cloglog = logit - cloglog,
         better_AIC    = case_when(logit_cloglog > 0 ~ "cloglog",
                                   logit_cloglog < 0 ~ "logit",
                                   TRUE              ~ "tied"))
link_diff_AIC %>% 
  count(better_AIC) # logit superior

ggplot(link_diff_AIC, aes(better_AIC, abs(logit_cloglog))) +
  geom_boxplot() # if logit better, larger difference

link_diff_AIC %>% 
  count(simtype, better_AIC) # logit superior

link_diff_AIC %>% 
  count(nat_haz, better_AIC) # logit superior

link_diff_AIC %>% 
  count(profile, better_AIC) # logit superior

link_diff_AIC %>% 
  count(resp_var_type, better_AIC) # logit superior

link_diff_AIC %>% 
  count(formula_nr, better_AIC) # logit superior

# keep logit models

glances_red <- glances %>% 
  filter(link == "logit")


# compare formulas --------------------------------------------------------
## AIC
ggplot(glances_red, aes(factor(formula_nr), AIC)) +
  geom_boxplot() # no apparent differences

formula_diff_AIC <- glances_red %>% 
  select(simtype:link, AIC) %>% 
  pivot_wider(names_from   = formula_nr,
              names_prefix = "f",
              values_from  = AIC) %>% 
  mutate(best_formula = case_when(f1 < f2 & f1 < f3 ~ "f1",
                                  f2 < f1 & f2 < f3 ~ "f2",
                                  f3 < f1 & f3 < f2 ~ "f3",
                                  TRUE              ~ "tie"))

formula_diff_AIC %>% 
  count(best_formula) # f1: 7; f2: 6; f3: 3

formula_diff_AIC %>% 
  count(simtype, best_formula) # LT: f1 best (5/0/3); ST: f2 best (2/6/0)

formula_diff_AIC %>% 
  count(nat_haz, best_formula) # A: f2 best (2/4/2); LED: f1 best (5/2/1)

formula_diff_AIC %>% 
  count(profile, best_formula) # MP: f1 & f2 best (4/4/0); IP: f1 & f3 best (3/2/3)

formula_diff_AIC %>% 
  count(resp_var_type, best_formula) # abs: f1 & f2 best (3/3/2); diff: f1 best (4/3/1)

# Option 1
# keep f1 for LT models and f2 for ST models
# -> easier comparison but not always the same model

# Option 2
# keep best model for each combination
# -> always the best model but comparison less straightforward




