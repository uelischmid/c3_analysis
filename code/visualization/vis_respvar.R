## visualization of response variable distributions
## 7.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
folder_in <- "data/raw/nais_analysis/"
folder_out <- "results/vis_respvar/"


# load data ---------------------------------------------------------------
LT <- read_rds(str_c(folder_in, "LT.rds")) %>% 
  mutate(simtype = "LT")
ST <- read_rds(str_c(folder_in, "ST.rds")) %>% 
  mutate(simtype = "ST")


# merge and pivot ---------------------------------------------------------
comb <- bind_rows(LT, ST) %>% 
  mutate(stratum = factor(stratum, levels = c("UM", "HM", "SA")),
         quality = factor(quality),
         init    = factor(init),
         mgm     = factor(mgm, levels = c("NOM", "STS", "GRS1", "GRS2", "CAB1", "CAB2", "SC")),
         nat_haz = factor(nat_haz),
         simtype = factor(simtype),
         q_site2 = case_when(q_site == 3 ~ "medium",
                             TRUE        ~ "good"),
         q_site2 = factor(q_site2, levels = c("medium", "good")),
         q_reg2  = case_when(q_reg == 1 ~ "bad",
                             TRUE       ~ "good"),
         q_reg2  = factor(q_reg2, levels = c("bad", "good"))) %>% 
  pivot_longer(cols = sha_y_MP_met:neg_dist_IP,
               names_to = "resp_var",
               values_to = "value")


# plot --------------------------------------------------------------------
# response variable distribution
ggplot(comb, aes(resp_var, value, color = simtype)) + 
  geom_boxplot() +
  labs(title = "Potential response variables",
       x     = "response variable") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

ggsave(str_c(folder_out, "0_resp_var_distr.png"))

# mean_i_total
data <- comb
respvar <- "mean_i_total"

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
  plots[[6]] <- ggplot(data_red, aes(mgm, value, color = simtype)) +
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


vars <- unique(comb$resp_var)

for (i in seq_along(vars)) {
  gg <- plot_explvar(respvar = vars[i],
                     data    = comb)
  ggsave(str_c(folder_out, vars[i], ".jpg"),
         plot = gg,
         scale = 1.5)
}
