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
  mutate(sha_i_MP_met_abs_logit = boot::logit(sha_i_MP_met_abs)) %>% 
  as.data.frame()

data_LT_A_UM <- data_LT_A %>% 
  filter(stratum == "UM")

data_LT_A_HM <- data_LT_A %>% 
  filter(stratum == "HM")

data_LT_A_SA <- data_LT_A %>% 
  filter(stratum == "SA")

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


# build model -------------------------------------------------------------
# LT A abs MP
m_LT_A_abs_MP <- gbm.step(data            = data_LT_A,
                          gbm.x           = c(3, 6, 8, 10, 11, 13),
                          gbm.y           = 15,
                          family          = "gaussian",
                          tree.complexity = 5,
                          learning.rate   = 0.05,
                          bag.fraction    = 0.75)

m_LT_A_abs_MP_UM <- gbm.step(data            = data_LT_A_UM,
                             gbm.x           = c(6, 8, 10, 11, 13),
                             gbm.y           = 15,
                             family          = "gaussian",
                             tree.complexity = 4,
                             learning.rate   = 0.01,
                             bag.fraction    = 0.75)

m_LT_A_abs_MP_HM <- gbm.step(data            = data_LT_A_HM,
                             gbm.x           = c(6, 8, 10, 11, 13),
                             gbm.y           = 15,
                             family          = "gaussian",
                             tree.complexity = 4,
                             learning.rate   = 0.01,
                             bag.fraction    = 0.75)

m_LT_A_abs_MP_SA <- gbm.step(data            = data_LT_A_SA,
                             gbm.x           = c(6, 8, 10, 11, 13),
                             gbm.y           = 15,
                             family          = "gaussian",
                             tree.complexity = 4,
                             learning.rate   = 0.01,
                             bag.fraction    = 0.75)


# logit-transformed response
m_LT_A_abs_MP_transf <- gbm.step(data            = data_LT_A,
                                 gbm.x           = c(3, 6, 8, 10, 11, 13),
                                 gbm.y           = 21,
                                 family          = "gaussian",
                                 tree.complexity = 5,
                                 learning.rate   = 0.05,
                                 bag.fraction    = 0.75)

m_LT_A_abs_MP_UM_transf <- gbm.step(data            = data_LT_A_UM,
                                    gbm.x           = c(6, 8, 10, 11, 13),
                                    gbm.y           = 21,
                                    family          = "gaussian",
                                    tree.complexity = 4,
                                    learning.rate   = 0.01,
                                    bag.fraction    = 0.75)

m_LT_A_abs_MP_HM_transf <- gbm.step(data            = data_LT_A_HM,
                                    gbm.x           = c(6, 8, 10, 11, 13),
                                    gbm.y           = 21,
                                    family          = "gaussian",
                                    tree.complexity = 4,
                                    learning.rate   = 0.01,
                                    bag.fraction    = 0.75)

m_LT_A_abs_MP_SA_transf <- gbm.step(data            = data_LT_A_SA,
                                    gbm.x           = c(6, 8, 10, 11, 13),
                                    gbm.y           = 21,
                                    family          = "gaussian",
                                    tree.complexity = 4,
                                    learning.rate   = 0.01,
                                    bag.fraction    = 0.75)

# summaries
summary(m_LT_A_abs_MP)
summary(m_LT_A_abs_MP_UM) # intensity moves 1 up
summary(m_LT_A_abs_MP_HM) # type moves 1 up, site moves 1 up
summary(m_LT_A_abs_MP_SA) # interval and intensity on top

# comparison to transformed
m_LT_A_abs_MP$cv.statistics$correlation.mean
m_LT_A_abs_MP_transf$cv.statistics$correlation.mean

m_LT_A_abs_MP_UM$cv.statistics$correlation.mean
m_LT_A_abs_MP_UM_transf$cv.statistics$correlation.mean

m_LT_A_abs_MP_HM$cv.statistics$correlation.mean
m_LT_A_abs_MP_HM_transf$cv.statistics$correlation.mean

m_LT_A_abs_MP_SA$cv.statistics$correlation.mean
m_LT_A_abs_MP_SA_transf$cv.statistics$correlation.mean


# plot all
gbm.plot(m_LT_A_abs_MP,
         n.plots     = 6,
         plot.layout = c(2, 3),
         write.title = FALSE)
gbm.plot(m_LT_A_abs_MP_UM,
         n.plots     = 5,
         plot.layout = c(2, 3),
         write.title = FALSE)
gbm.plot(m_LT_A_abs_MP_HM,
         n.plots     = 5,
         plot.layout = c(2, 3),
         write.title = FALSE)
gbm.plot(m_LT_A_abs_MP_SA,
         n.plots     = 5,
         plot.layout = c(2, 3),
         write.title = FALSE)


gbm.plot(m_LT_A_abs_MP_transf,
         n.plots     = 6,
         plot.layout = c(2, 3),
         write.title = FALSE)
source("code/visualization/gbm.plot_us.R")
gbm.plot_us_invlogit(m_LT_A_abs_MP_transf,
                     n.plots     = 6,
                     plot.layout = c(2, 3),
                     write.title = FALSE)


gbm.plot(m_LT_A_abs_MP_UM_transf,
         n.plots     = 5,
         plot.layout = c(2, 3),
         write.title = FALSE)
gbm.plot(m_LT_A_abs_MP_HM_transf,
         n.plots     = 5,
         plot.layout = c(2, 3),
         write.title = FALSE)
gbm.plot(m_LT_A_abs_MP_SA_transf,
         n.plots     = 5,
         plot.layout = c(2, 3),
         write.title = FALSE)



m_int <- gbm.interactions(m_LT_A_abs_MP)
m_int$rank.list
m_int$interactions

m_int_UM <- gbm.interactions(m_LT_A_abs_MP_UM)
m_int_UM$rank.list

m_int_HM <- gbm.interactions(m_LT_A_abs_MP_HM)
m_int_HM$rank.list

m_int_SA <- gbm.interactions(m_LT_A_abs_MP_SA)
m_int_SA$rank.list

gbm.perspec(m_LT_A_abs_MP,
            x = 1,
            y = 5,
            z.range = c(0.1, 0.7))
gbm.perspec(m_LT_A_abs_MP,
            x = 5,
            y = 6,
            z.range = c(0, 0.3))

gbm.perspec(m_LT_A_abs_MP_UM,
            x = 5,
            y = 4,
            z.range = c(0, 0.6))

gbm.perspec(m_LT_A_abs_MP_UM,
            x = 3,
            y = 4,
            z.range = c(0, 0.6))

gbm.perspec(m_LT_A_abs_MP_HM,
            x = 5,
            y = 4,
            z.range = c(0, 0.3))

gbm.perspec(m_LT_A_abs_MP_SA,
            x = 5,
            y = 4,
            z.range = c(0.4, 0.7))

# for (i in seq_along(data_LT_A)) {
#   cat(i, ": ", colnames(data_LT_A)[i], "\n")
# }
