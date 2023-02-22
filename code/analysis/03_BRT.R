## analyze data with boosted regression trees
## 22.2.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)
library(dismo)

folder_in <- "data/processed/analysis/"
folder_out <- "results/brt_models/"


# load data ---------------------------------------------------------------
data_analysis <- read_rds(str_c(folder_in, "analysis_data_transf.rds"))

# LT A
data_LT_A <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "A") %>% 
  as.data.frame()

# data_LT_A_UM <- data_LT_A %>% 
#   filter(stratum == "UM")
# 
# data_LT_A_HM <- data_LT_A %>% 
#   filter(stratum == "HM")
# 
# data_LT_A_SA <- data_LT_A %>% 
#   filter(stratum == "SA")

# LT LED
data_LT_LED <- data_analysis %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "LED") %>% 
  as.data.frame()

# ST A
data_ST_A <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "A") %>% 
  as.data.frame()

# ST LED
data_ST_LED <- data_analysis %>% 
  filter(simtype == "ST") %>% 
  filter(nat_haz == "LED") %>% 
  as.data.frame()


# build models ------------------------------------------------------------


# for (i in seq_along(data_LT_A)) {
#   cat(i, ": ", colnames(data_LT_A)[i], "\n")
# }
