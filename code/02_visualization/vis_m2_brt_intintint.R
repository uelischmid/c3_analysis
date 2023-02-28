## visualize boosted regression trees
## intensity - interval - interaction
## 28.2.23, us


# setup -------------------------------------------------------------------
library(tidyverse)
library(gbm)

folder_in <- "data/processed/brt_models_v2/"
folder_out <- "results/vis_models_v2/brt/vis_brt_red_intintint/"


# load data ---------------------------------------------------------------
models_LT <- read_rds(str_c(folder_in, "models_red_LT.rds"))
model_combinations_LT <- read_rds(str_c(folder_in, "model_combinations_red_LT.rds")) %>% 
  rowid_to_column()

models_ST <- read_rds(str_c(folder_in, "models_red_ST.rds"))
model_combinations_ST <- read_rds(str_c(folder_in, "model_combinations_red_ST.rds")) %>% 
  rowid_to_column()


# function ----------------------------------------------------------------
gbm_perspec_int3 <- function(gbm.object,
                             plot.title,
                             x          = 4,
                             x.label    = "Intensity (%)",
                             y          = 3,
                             y.label    = "Interval (y)",
                             z.label    = "Response",
                             pred.means = NULL,
                             theta      = 55,
                             phi        = 40,
                             ticktype   = "detailed") {
  
  gbm.call <- gbm.object$gbm.call
  gbm.x <- gbm.call$gbm.x
  n.preds <- length(gbm.x)
  pred.names <- gbm.call$predictor.names

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
  
  z.range <- c(floor(min(prediction) * 10) / 10,
               ceiling(max(prediction) * 10) / 10)
  
  pred.matrix <- matrix(prediction, ncol = 50, nrow = 50)
  
  persp(x        = x.var,
        y        = y.var,
        z        = pred.matrix,
        zlim     = z.range, 
        xlab     = x.label,
        ylab     = y.label,
        zlab     = z.label, 
        theta    = theta,
        phi      = phi,
        r        = sqrt(10),
        d        = 3, 
        ticktype = ticktype,
        mgp      = c(4, 1, 0),
        main     = plot.title)
}


# plot --------------------------------------------------------------------
# LT
for (i in 1:nrow(model_combinations_LT)) {
  cat("\r", i, "/", nrow(model_combinations_LT))
  vals <- model_combinations_LT[i,]
  plot_title <- str_c(pull(vals, simtype),
                      pull(vals, nat_haz),
                      pull(vals, profile),
                      pull(vals, resp_var_type),
                      pull(vals, stratum),
                      sep = " ")
  file_name <- str_c(str_replace_all(plot_title, " ", "_"), ".jpg")
  
  jpeg(filename  = str_c(folder_out, file_name),
       width     = 1000,
       height    = 1000,
       quality   = 150,
       pointsize = 20)
  gbm_perspec_int3(gbm.object = models_LT[[i]],
                   plot.title = plot_title)
  dev.off()
}


# ST
for (i in 1:nrow(model_combinations_ST)) {
  cat("\r", i, "/", nrow(model_combinations_ST))
  vals <- model_combinations_ST[i,]
  plot_title <- str_c(pull(vals, simtype),
                      pull(vals, nat_haz),
                      pull(vals, profile),
                      pull(vals, resp_var_type),
                      pull(vals, stratum),
                      pull(vals, init),
                      sep = " ")
  file_name <- str_c(str_replace_all(plot_title, " ", "_"), ".jpg")
  
  jpeg(filename  = str_c(folder_out, file_name),
       width     = 1000,
       height    = 1000,
       quality   = 150,
       pointsize = 20)
  gbm_perspec_int3(gbm.object = models_ST[[i]],
                   plot.title = plot_title)
  dev.off()
}
