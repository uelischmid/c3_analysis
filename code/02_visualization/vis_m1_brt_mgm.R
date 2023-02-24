## visualize boosted regression trees
## management effects
## 24.2.23, us

# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)
library(dismo)
library(cowplot)
library(gridGraphics)

folder_in <- "data/processed/brt_models_v1/"
folder_out <- "results/vis_models_v1/brt/vis_brt_mgm/"


# load data ---------------------------------------------------------------
# models_full <- read_rds(str_c(folder_in, "models_full.rds"))
# model_combinations_full <- read_rds(str_c(folder_in, "model_combinations_full.rds")) %>% 
#   rowid_to_column()
models_stratum <- read_rds(str_c(folder_in, "models_stratum.rds"))
model_combinations_stratum <- read_rds(str_c(folder_in, "model_combinations_stratum.rds")) %>% 
  rowid_to_column()


# functions ---------------------------------------------------------------
get_interaction_zminmax <- function(gbm.object, x, y, pred.means = NULL) {
  gbm.call <- gbm.object$gbm.call
  gbm.x <- gbm.call$gbm.x
  n.preds <- length(gbm.x)
  gbm.y <- gbm.call$gbm.y
  pred.names <- gbm.call$predictor.names
  family = gbm.call$family
  have.factor <- FALSE
  
  x.name <- gbm.call$predictor.names[x]
  y.name <- gbm.call$predictor.names[y]
  
  data <- gbm.call$dataframe[, gbm.x, drop = FALSE]
  n.trees <- gbm.call$best.trees
  
  x.var <- seq(min(data[, x], na.rm = T), max(data[, x], na.rm = T), length = 50)
  y.var <- seq(min(data[, y], na.rm = T), max(data[, y], na.rm = T), length = 50)
  
  pred.frame <- expand.grid(list(x.var, y.var))
  names(pred.frame) <- c(x.name, y.name)
  pred.rows <- nrow(pred.frame)
  
  j <- 3
  for (i in 1:n.preds) {
    if (i != x & i != y) {
      if (is.vector(data[, i])) {
        m <- match(pred.names[i], names(pred.means))
        if (is.na(m)) {
          pred.frame[, j] <- mean(data[, i], na.rm = T)
        } else {
          pred.frame[, j] <- pred.means[m]
        }
      }
      
      if (is.factor(data[, i])) {
        m <- match(pred.names[i], names(pred.means))
        temp.table <- table(data[, i])
        if (is.na(m)) {
          pred.frame[, j] <- rep(names(temp.table)[2], pred.rows)
        } else {
          pred.frame[, j] <- pred.means[m]
        }
        pred.frame[, j] <- factor(pred.frame[, j], levels = names(temp.table))
      }
      
      names(pred.frame)[j] <- pred.names[i]
      j <- j + 1
    }
  }
  prediction <- gbm::predict.gbm(gbm.object, pred.frame, n.trees = n.trees, 
                                 type = "response")
  
  minmax <- c(min(prediction), max(prediction))
  return(minmax)
}


# plot --------------------------------------------------------------------
for (i in 1:nrow(model_combinations_stratum)) {
  vals <- model_combinations_stratum[i, ]
  
  m_used <- models_stratum[[pull(vals, rowid)]]
  vno_type <- match("mgm_type", m_used$var.names)
  vno_intensity <- match("mgm_intensity", m_used$var.names)
  vno_interval <- match("mgm_interval", m_used$var.names)
  
  p_type <- ~gbm.plot(gbm.object   = m_used,
                      variable.no  = vno_type,
                      plot.layout  = c(1, 1),
                      x.label      = "Type",
                      y.label      = "Response",
                      write.title  = FALSE,
                      show.contrib = FALSE)
  
  zrange_pred <- get_interaction_zminmax(gbm.object = m_used,
                                         x          = vno_intensity,
                                         y          = vno_interval)
  z_range_used <- c(floor(zrange_pred[1] * 10) / 10,
                    ceiling(zrange_pred[2] * 10) / 10)
  
  p_intint <- ~gbm.perspec(gbm.object = m_used,
                           x          = vno_intensity,
                           x.label    = "Intensity (%)",
                           y          = vno_interval,
                           y.label    = "Interval (y)",
                           z.label    = "Response",
                           z.range    = z_range_used)
  
  p_both <- plot_grid(p_type, p_intint,
                      ncol = 2)
  ggsave2(str_c(folder_out,
                pull(vals, simtype), "_",
                pull(vals, nat_haz), "_",
                pull(vals, profile), "_",
                pull(vals, resp_var_type), "_",
                pull(vals, stratum), ".jpg"),
          p_both,
          height = 15,
          width = 30,
          units = "cm")
}
