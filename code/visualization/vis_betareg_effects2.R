## visualize effects of betareg-models
## 20.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)
# library(effects)
library(ggeffects)

folder_in_d <- "data/processed/analysis/"
folder_in_m <- "results/betareg_models/"
folder_out <- "results/vis_betareg_effects/"


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
  filter(nat_haz == "A")

data_ST_LED <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "LED")


# load models -------------------------------------------------------------
models <- read_rds(str_c(folder_in_m, "models.rds"))
model_combinations <- read_rds(str_c(folder_in_m, "model_combinations.rds")) %>% 
  rowid_to_column()



# visualize LT abs --------------------------------------------------------
mc_used <- model_combinations %>% 
  filter(link == "logit") %>% 
  filter(simtype == "ST") %>% 
  filter(resp_var_type == "abs") %>% 
  filter(formula_nr == 3)

summary(models[[53]])
plot(effects::allEffects(models[[53]]))

# plot stratum
effs <- vector(mode = "list", nrow(mc_used))
effs_df <- effs
for (i in 1:nrow(mc_used)) {
  vals <- mc_used[i,]
  effs[[i]] <- ggeffect(model = models[[pull(vals, rowid)]],
                        terms = "stratum")
  
  effs_df[[i]] <- effs[[i]] %>% 
    as_tibble() %>% 
    mutate(nat_haz = pull(vals, nat_haz),
           profile = pull(vals, profile))
}
effs_df <- bind_rows(effs_df)

ggplot(effs_df, aes(x, predicted, color = nat_haz)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0.25,
                position = position_dodge(width = 0.2)) +
  facet_grid(~profile) +
  labs(x = get_x_title(effs[[1]]),
       y = get_y_title(effs[[1]]))


# plot qualities
effs <- vector(mode = "list", nrow(mc_used))
effs_df <- effs
for (i in 1:nrow(mc_used)) {
  vals <- mc_used[i,]
  effs[[i]] <- ggeffect(model = models[[pull(vals, rowid)]],
                        terms = c("q_site2", "q_reg2"))
  
  effs_df[[i]] <- effs[[i]] %>% 
    as_tibble() %>% 
    mutate(nat_haz = pull(vals, nat_haz),
           profile = pull(vals, profile))
}
effs_df <- bind_rows(effs_df) %>% 
  rename(q_reg2 = group)

ggplot(effs_df, aes(x, predicted, color = nat_haz, shape = q_reg2, linetype = q_reg2)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0.25,
                position = position_dodge(width = 0.2)) +
  facet_wrap(~profile) +
  labs(x = get_x_title(effs[[1]]),
       y = get_y_title(effs[[1]]))


# plot mgm type
effs <- vector(mode = "list", nrow(mc_used))
effs_df <- effs
for (i in 1:nrow(mc_used)) {
  vals <- mc_used[i,]
  effs[[i]] <- ggeffect(model = models[[pull(vals, rowid)]],
                        terms = "mgm_type")
  
  effs_df[[i]] <- effs[[i]] %>% 
    as_tibble() %>% 
    mutate(nat_haz = pull(vals, nat_haz),
           profile = pull(vals, profile))
}
effs_df <- bind_rows(effs_df)

ggplot(effs_df, aes(x, predicted, color = nat_haz)) +
  geom_point(position = position_dodge(width = 0.2)) +
  geom_errorbar(aes(ymin = conf.low,
                    ymax = conf.high),
                width = 0.25,
                position = position_dodge(width = 0.2)) +
  facet_grid(~profile) +
  labs(x = get_x_title(effs[[1]]),
       y = get_y_title(effs[[1]]))


# plot mgm interval
effs <- vector(mode = "list", nrow(mc_used))
effs_df <- effs
for (i in 1:nrow(mc_used)) {
  vals <- mc_used[i,]
  effs[[i]] <- ggeffect(model = models[[pull(vals, rowid)]],
                        terms = "mgm_interval [10:40]") 
  
  effs_df[[i]] <- effs[[i]] %>% 
    as_tibble() %>% 
    mutate(nat_haz = pull(vals, nat_haz),
           profile = pull(vals, profile))
}

effs_df <- bind_rows(effs_df)

ggplot(effs_df, aes(x, predicted, color = nat_haz)) +
  geom_ribbon(aes(x, predicted,
                  ymin = conf.low, ymax = conf.high,
                  fill = nat_haz),
              alpha = 0.2, color = NA,
              inherit.aes = FALSE) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~profile) +
  labs(x = get_x_title(effs[[1]]),
       y = get_y_title(effs[[1]]))


# plot mgm_intensity
effs <- vector(mode = "list", nrow(mc_used))
effs_df <- effs
for (i in 1:nrow(mc_used)) {
  vals <- mc_used[i,]
  effs[[i]] <- ggeffect(model = models[[pull(vals, rowid)]],
                        terms = "mgm_intensity [10:40]") 
  
  effs_df[[i]] <- effs[[i]] %>% 
    as_tibble() %>% 
    mutate(nat_haz = pull(vals, nat_haz),
           profile = pull(vals, profile))
}

effs_df <- bind_rows(effs_df)

ggplot(effs_df, aes(x, predicted, color = nat_haz)) +
  geom_ribbon(aes(x, predicted,
                  ymin = conf.low, ymax = conf.high,
                  fill = nat_haz),
              alpha = 0.2, color = NA,
              inherit.aes = FALSE) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~profile) +
  labs(x = get_x_title(effs[[1]]),
       y = get_y_title(effs[[1]]))

# plot int-int-interaction
effs <- vector(mode = "list", nrow(mc_used))
effs_df <- effs
for (i in 1:nrow(mc_used)) {
  vals <- mc_used[i,]
  effs[[i]] <- ggeffect(model = models[[pull(vals, rowid)]],
                        terms = c("mgm_intensity [10:40]",
                                  "mgm_interval [10, 20, 30, 40]")) 
  
  effs_df[[i]] <- effs[[i]] %>% 
    as_tibble() %>% 
    mutate(nat_haz = pull(vals, nat_haz),
           profile = pull(vals, profile))
}

effs_df <- bind_rows(effs_df) %>% 
  rename(mgm_interval = group)


effs_df %>% 
  filter(profile == "MP") %>% 
  ggplot(aes(x, predicted, color = nat_haz)) +
  geom_ribbon(aes(x, predicted,
                  ymin = conf.low, ymax = conf.high,
                  fill = nat_haz),
              alpha = 0.2, color = NA,
              inherit.aes = FALSE) +
  geom_line() +
  facet_wrap(~mgm_interval,
             labeller = "label_both") +
  labs(title = "MP",
       x = get_x_title(effs[[1]]),
       y = get_y_title(effs[[1]]))

effs_df %>% 
  filter(profile == "IP") %>% 
  ggplot(aes(x, predicted, color = nat_haz)) +
  geom_ribbon(aes(x, predicted,
                  ymin = conf.low, ymax = conf.high,
                  fill = nat_haz),
              alpha = 0.2, color = NA,
              inherit.aes = FALSE) +
  geom_line() +
  facet_wrap(~mgm_interval,
             labeller = "label_both") +
  labs(title = "IP",
       x = get_x_title(effs[[1]]),
       y = get_y_title(effs[[1]]))
