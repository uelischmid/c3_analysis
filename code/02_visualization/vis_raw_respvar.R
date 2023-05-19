## visualization of response variable distributions
## 9.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(cowplot)

folder_in <- "data/processed/nais_analysis_data/"
folder_out <- "results/vis_raw_data/vis_respvar/"


# load data ---------------------------------------------------------------
comb <- read_rds(str_c(folder_in, "analysis_data_transf.rds"))
respvar_names <- colnames(comb)[15:18]
comb <- comb %>% 
  select(simtype:sha_i_IP_met_diff) %>% 
  pivot_longer(cols      = sha_i_MP_met_abs:sha_i_IP_met_diff,
               names_to  = "resp_var",
               values_to = "value") %>% 
  mutate(resp_var = factor(resp_var, levels = respvar_names))


# plot response variable distributions ------------------------------------
# all response variables
ggplot(comb, aes(resp_var, value, color = simtype)) + 
  geom_boxplot() +
  labs(title = "Potential response variables",
       x     = "response variable") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(str_c(folder_out, "0_resp_var_distr_all.png"))

# selected response variables (boxplots)
comb_sel <- comb %>% 
  filter(str_detect(resp_var, "sha_i_"))

ggplot(comb_sel, aes(resp_var, value, color = simtype)) + 
  geom_boxplot() +
  labs(title = "Selected response variables",
       x     = "response variable") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(str_c(folder_out, "0_resp_var_distr_sel_1.png"))

# selected response variables (density)
p1 <- comb_sel %>% 
  filter(str_detect(resp_var, "abs")) %>% 
  ggplot(aes(value)) +
  geom_density() +
  geom_vline(aes(xintercept = min),
             data = comb_sel %>% 
               filter(str_detect(resp_var, "abs")) %>% 
               group_by(simtype, resp_var) %>% 
               summarise(min = min(value)),
             lty = 2) +
  geom_vline(aes(xintercept = max),
             data = comb_sel %>% 
               filter(str_detect(resp_var, "abs")) %>% 
               group_by(simtype, resp_var) %>% 
               summarise(max = max(value)),
             lty = 2) +
  facet_grid(rows = vars(simtype),
             cols = vars(resp_var)) +
  xlim(0, 1) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

p2 <- comb_sel %>% 
  filter(str_detect(resp_var, "diff")) %>% 
  ggplot(aes(value)) +
  geom_density() +
  geom_vline(aes(xintercept = min),
             data = comb_sel %>% 
               filter(str_detect(resp_var, "diff")) %>% 
               group_by(simtype, resp_var) %>% 
               summarise(min = min(value)),
             lty = 2) +
  geom_vline(aes(xintercept = max),
             data = comb_sel %>% 
               filter(str_detect(resp_var, "diff")) %>% 
               group_by(simtype, resp_var) %>% 
               summarise(max = max(value)),
             lty = 2) +
  facet_grid(rows = vars(simtype),
             cols = vars(resp_var)) +
  xlim(-1, 1) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

plot_grid(p1, p2,
                ncol = 1)
ggsave(str_c(folder_out, "0_resp_var_distr_sel_2.png"))




# plot explanatory variables ----------------------------------------------
plot_explvar <- function(respvar, data) {
  plots <- vector(mode = "list", length = 10)
  
  data_red <- data %>%
    filter(resp_var == respvar)
  
  # stratum
  plots[[1]] <- ggplot(data_red, aes(stratum, value, color = simtype)) +
    geom_boxplot() +
    labs(title = "stratum",
         y     = respvar) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # quality
  plots[[2]] <- ggplot(data_red, aes(quality, value, color = simtype)) +
    geom_boxplot() +
    labs(title = "combined quality",
         y     = respvar) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # site quality
  plots[[3]] <- ggplot(data_red, aes(q_site2, value, color = simtype)) +
    geom_boxplot() +
    labs(title = "site quality",
         y     = respvar) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # regeneration quality
  plots[[4]] <- ggplot(data_red, aes(q_reg2, value, color = simtype)) +
    geom_boxplot() +
    labs(title = "regeneration quality",
         y     = respvar) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # init
  plots[[5]] <- ggplot(data_red, aes(init, value, color = simtype)) +
    geom_boxplot() +
    labs(title = "init",
         y     = respvar) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # mgm type
  plots[[6]] <- ggplot(data_red, aes(mgm_type, value, color = simtype)) +
    geom_boxplot() +
    labs(title = "management type",
         y     = respvar) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # mgm interval
  plots[[7]] <- ggplot(data_red, aes(factor(mgm_interval), value, color = simtype)) +
    geom_boxplot() +
    labs(title = "management interval",
         x     = "management interval",
         y     = respvar) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # mgm intensity
  plots[[8]] <- ggplot(data_red, aes(factor(mgm_intensity), value, color = simtype)) +
    geom_boxplot() +
    labs(title = "management intensity",
         x     = "management intensity",
         y     = respvar) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # natural hazard
  plots[[9]] <- ggplot(data_red, aes(nat_haz, value, color = simtype)) +
    geom_boxplot() +
    labs(title = "natural hazard",
         y     = respvar) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          legend.position = "none")
  
  # legend
  plots[[10]] <- get_legend(ggplot(data_red, aes(nat_haz, value, color = simtype)) +
               geom_boxplot())
  
  # assemble
  gg <- plot_grid(plotlist = plots)
  
  return(gg)
}

vars <- as.character(unique(comb_sel$resp_var))

for (i in seq_along(vars)) {
  gg <- plot_explvar(respvar = vars[i],
                     data    = comb_sel)
  ggsave(str_c(folder_out, vars[i], ".jpg"),
         plot = gg,
         scale = 1.5)
}
