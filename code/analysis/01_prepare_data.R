## prepare data
## calculate alternative response variables
## stanardise explanatory variables
## 9.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)

folder_in <- "data/raw/nais_analysis/"
folder_out <- "data/processed/analysis/"


# load data ---------------------------------------------------------------
LT_orig <- read_rds(str_c(folder_in, "LT.rds")) %>% 
  mutate(simtype = "LT")
ST_orig <- read_rds(str_c(folder_in, "ST.rds")) %>% 
  mutate(simtype = "ST")

comb_orig <- bind_rows(LT_orig, ST_orig)


# calculate alternative response variable (diff NOM) ----------------------
comb_mgm <- comb_orig %>% 
  filter(mgm != "NOM") %>% 
  rename_with(~ str_c(.x, "_abs"), c(sha_y_MP_met:neg_dist_IP))
comb_nom <- comb_orig %>% 
  filter(mgm == "NOM") %>% 
  rename_with(~ str_c(.x, "_NOM"), c(sha_y_MP_met:neg_dist_IP)) %>% 
  select(-c(mgm, mgm_interval, mgm_intensity))

comb <- left_join(comb_mgm, comb_nom, by = c("stratum",
                                             "quality", "q_site", "q_reg",
                                             "init", "nat_haz", "simtype")) %>% 
  mutate(sha_y_MP_met_diff = sha_y_MP_met_abs - sha_y_MP_met_NOM,
         sha_y_IP_met_diff = sha_y_IP_met_abs - sha_y_IP_met_NOM,
         sha_i_MP_met_diff = sha_i_MP_met_abs - sha_i_MP_met_NOM,
         sha_i_IP_met_diff = sha_i_IP_met_abs - sha_i_IP_met_NOM,
         mean_i_total_diff = mean_i_total_abs - mean_i_total_NOM,
         neg_dist_MP_diff  = neg_dist_MP_abs - neg_dist_MP_NOM,
         neg_dist_IP_diff  = neg_dist_IP_abs - neg_dist_IP_NOM) %>% 
  select(simtype, nat_haz, stratum:mgm_intensity,
         sha_i_MP_met_abs, sha_i_MP_met_diff,
         sha_i_IP_met_abs, sha_i_IP_met_diff,
         neg_dist_MP_abs, neg_dist_MP_diff,
         neg_dist_IP_abs, neg_dist_IP_diff,
         mean_i_total_abs, mean_i_total_diff,
         sha_y_MP_met_abs, sha_y_MP_met_diff,
         sha_y_IP_met_abs, sha_y_IP_met_diff)




# convert response variables ----------------------------------------------
# standardize mgm interval & intensity
comb <- comb %>% 
  mutate(mgm_interval_s  = as.vector(scale(mgm_interval)),
         mgm_intensity_s = as.vector(scale(mgm_intensity)))


# factors
comb <- comb %>% 
  mutate(across(simtype:mgm, ~ factor(.x))) %>% 
  mutate(stratum = fct_relevel(stratum, "UM", "HM", "SA"),
         mgm     = fct_relevel(mgm, "STS", "GRS1", "GRS2", "CAB1", "CAB2", "SC"),
         q_site2 = q_site,
         q_site2 = fct_recode(q_site2,
                              medium = "3",
                              good   = "5"),
         q_site2 = fct_relevel(q_site2, "good", "medium"),
         q_reg2  = q_reg,
         q_reg2  = fct_recode(q_reg2,
                              hindered = "1",
                              normal   = "3",
                              normal   = "5"),
         q_reg2  = fct_relevel(q_reg2, "normal", "hindered"))

# reorder
comb <- comb %>% 
  select(simtype, nat_haz, stratum,
         quality, q_site, q_site2, q_reg, q_reg2,
         init,
         mgm_type = mgm, mgm_interval, mgm_interval_s, mgm_intensity, mgm_intensity_s,
         sha_i_MP_met_abs:sha_y_IP_met_diff)

# save --------------------------------------------------------------------
write_rds(comb,
          str_c(folder_out, "analysis_data_transf.rds"))


