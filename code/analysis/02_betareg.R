## analyize data with beta regression
## 9.2.23, us


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


# build models for absolute values ----------------------------------------
abs_model_info <- vector(mode = "list", length = 8)

## LT A MP
br_LT_A_MP_abs <- betareg(sha_i_MP_met_abs ~ stratum +
                            q_site2 + q_reg2 +
                            mgm_type +
                            mgm_interval_s + mgm_intensity_s +
                            I(mgm_interval_s^2) + I(mgm_intensity_s^2) +
                            q_site2 * q_reg2 +
                            mgm_interval_s * mgm_intensity_s,
                          data = data_LT_A)
summary(br_LT_A_MP_abs)
abs_model_info[[1]] <- glance(br_LT_A_MP_abs) %>% 
  mutate(model = "br_LT_A_MP_abs",
         simtype = "LT",
         nat_haz = "A",
         profile = "MP")

# LT A IP
br_LT_A_IP_abs <- betareg(sha_i_IP_met_abs ~ stratum +
                            q_site2 + q_reg2 +
                            mgm_type +
                            mgm_interval_s + mgm_intensity_s +
                            I(mgm_interval_s^2) + I(mgm_intensity_s^2) +
                            q_site2 * q_reg2 +
                            mgm_interval_s * mgm_intensity_s,
                          data = data_LT_A)
summary(br_LT_A_IP_abs)
abs_model_info[[2]] <- glance(br_LT_A_IP_abs) %>% 
  mutate(model = "br_LT_A_IP_abs",
         simtype = "LT",
         nat_haz = "A",
         profile = "IP")

# LT LED MP
br_LT_LED_MP_abs <- betareg(sha_i_MP_met_abs ~ stratum +
                            q_site2 + q_reg2 +
                            mgm_type +
                            mgm_interval_s + mgm_intensity_s +
                            I(mgm_interval_s^2) + I(mgm_intensity_s^2) +
                            q_site2 * q_reg2 +
                            mgm_interval_s * mgm_intensity_s,
                          data = data_LT_LED)
summary(br_LT_LED_MP_abs)
abs_model_info[[3]] <- glance(br_LT_LED_MP_abs) %>% 
  mutate(model = "br_LT_LED_MP_abs",
         simtype = "LT",
         nat_haz = "LED",
         profile = "MP")

# LT LED IP
br_LT_LED_IP_abs <- betareg(sha_i_IP_met_abs ~ stratum +
                            q_site2 + q_reg2 +
                            mgm_type +
                            mgm_interval_s + mgm_intensity_s +
                            I(mgm_interval_s^2) + I(mgm_intensity_s^2) +
                            q_site2 * q_reg2 +
                            mgm_interval_s * mgm_intensity_s,
                          data = data_LT_LED)
summary(br_LT_LED_IP_abs)
abs_model_info[[4]] <- glance(br_LT_LED_IP_abs) %>% 
  mutate(model = "br_LT_LED_IP_abs",
         simtype = "LT",
         nat_haz = "LED",
         profile = "IP")

## ST A MP
br_ST_A_MP_abs <- betareg(sha_i_MP_met_abs ~ stratum +
                            q_site2 + q_reg2 +
                            init +
                            mgm_type +
                            mgm_interval_s + mgm_intensity_s +
                            I(mgm_interval_s^2) + I(mgm_intensity_s^2) +
                            q_site2 * q_reg2 +
                            mgm_interval_s * mgm_intensity_s,
                          data = data_ST_A)
summary(br_ST_A_MP_abs)
abs_model_info[[5]] <- glance(br_ST_A_MP_abs) %>% 
  mutate(model = "br_ST_A_MP_abs",
         simtype = "ST",
         nat_haz = "A",
         profile = "MP")

# ST A IP
br_ST_A_IP_abs <- betareg(sha_i_IP_met_abs ~ stratum +
                            q_site2 + q_reg2 +
                            init +
                            mgm_type +
                            mgm_interval_s + mgm_intensity_s +
                            I(mgm_interval_s^2) + I(mgm_intensity_s^2) +
                            q_site2 * q_reg2 +
                            mgm_interval_s * mgm_intensity_s,
                          data = data_ST_A)
summary(br_ST_A_IP_abs)
abs_model_info[[6]] <- glance(br_ST_A_IP_abs) %>% 
  mutate(model = "br_ST_A_IP_abs",
         simtype = "ST",
         nat_haz = "A",
         profile = "IP")

# ST LED MP
br_ST_LED_MP_abs <- betareg(sha_i_MP_met_abs ~ stratum +
                              q_site2 + q_reg2 +
                              init +
                              mgm_type +
                              mgm_interval_s + mgm_intensity_s +
                              I(mgm_interval_s^2) + I(mgm_intensity_s^2) +
                              q_site2 * q_reg2 +
                              mgm_interval_s * mgm_intensity_s,
                            data = data_ST_LED)
summary(br_ST_LED_MP_abs)
abs_model_info[[7]] <- glance(br_ST_LED_MP_abs) %>% 
  mutate(model = "br_ST_LED_MP_abs",
         simtype = "ST",
         nat_haz = "LED",
         profile = "MP")

# ST LED IP
br_ST_LED_IP_abs <- betareg(sha_i_IP_met_abs ~ stratum +
                              q_site2 + q_reg2 +
                              init +
                              mgm_type +
                              mgm_interval_s + mgm_intensity_s +
                              I(mgm_interval_s^2) + I(mgm_intensity_s^2) +
                              q_site2 * q_reg2 +
                              mgm_interval_s * mgm_intensity_s,
                            data = data_ST_LED)
summary(br_ST_LED_IP_abs)
abs_model_info[[8]] <- glance(br_ST_LED_IP_abs) %>% 
  mutate(model = "br_ST_LED_IP_abs",
         simtype = "ST",
         nat_haz = "LED",
         profile = "IP")


# overview ----------------------------------------------------------------
abs_model_info <- bind_rows(abs_model_info) %>% 
  select(model:profile, pseudo.r.squared:nobs)


# visualize effects -------------------------------------------------------
plot(allEffects(br_LT_A_MP_abs))

model_names <- ls()[str_detect(ls(), "br_")]

for (i in seq_along(model_names)) {
  jpeg(filename = str_c(folder_out, model_names[i], ".jpg"),
       width = 1000,
       height = 700) 
  plot(allEffects(mod = get(model_names[i])))
  dev.off()
}
