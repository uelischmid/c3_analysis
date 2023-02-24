# gbm.object = mymodel
# x = 5
# y = 1
# pred.means = NULL
# x.label = NULL
# x.range = NULL
# y.label = NULL
# z.label = "fitted value"
# y.range = NULL
# z.range = c(0, 1)
# leg.coords = NULL
# ticktype = "detailed"
# theta = 55
# phi = 40
# smooth = "none"
# mask = FALSE
# perspective = TRUE

gbm.perspec_us <- function (gbm.object, x = 1, y = 2, pred.means = NULL, x.label = NULL, 
          x.range = NULL, y.label = NULL, z.label = "fitted value", 
          y.range = NULL, z.range = NULL, leg.coords = NULL, ticktype = "detailed", 
          theta = 55, phi = 40, smooth = "none", mask = FALSE, 
          perspective = TRUE, ...) {
  if (!requireNamespace("gbm")) {
    stop("you need to install the gbm package to use this function")
  }
  requireNamespace("splines")
  gbm.call <- gbm.object$gbm.call
  gbm.x <- gbm.call$gbm.x
  n.preds <- length(gbm.x)
  gbm.y <- gbm.call$gbm.y
  pred.names <- gbm.call$predictor.names
  family = gbm.call$family
  have.factor <- FALSE
  
  x.name <- gbm.call$predictor.names[x]
  if (is.null(x.label)) {
    x.label <- gbm.call$predictor.names[x]
  }
  
  y.name <- gbm.call$predictor.names[y]
  if (is.null(y.label)) {
    y.label <- gbm.call$predictor.names[y]
  }
  
  data <- gbm.call$dataframe[, gbm.x, drop = FALSE]
  n.trees <- gbm.call$best.trees
  
  if (is.vector(data[, x])) {
    if (is.null(x.range)) {
      x.var <- seq(min(data[, x], na.rm = T), max(data[, x], na.rm = T), length = 50)
    } else {
      x.var <- seq(x.range[1], x.range[2], length = 50)
    }
  } else {
    x.var <- names(table(data[, x]))
    have.factor <- TRUE
  }
  
  if (is.vector(data[, y])) {
    if (is.null(y.range)) {
      y.var <- seq(min(data[, y], na.rm = T), max(data[, 
                                                       y], na.rm = T), length = 50)
    } else {
      y.var <- seq(y.range[1], y.range[2], length = 50)
    }
  } else {
    y.var <- names(table(data[, y]))
    if (have.factor) {
      stop("at least one marginal predictor must be a vector!")
    } else {
      have.factor <- TRUE
    }
  }
  
  pred.frame <- expand.grid(list(x.var, y.var))
  names(pred.frame) <- c(x.name, y.name)
  pred.rows <- nrow(pred.frame)
  if (have.factor) {
    if (is.factor(pred.frame[, 2])) {
      pred.frame <- pred.frame[, c(2, 1)]
      x.var <- y.var
      y.label <- x.label # new
    }
  }
  
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
          pred.frame[, j] <- rep(names(temp.table)[2], 
                                 pred.rows)
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
  if (smooth == "model") {
    pred.glm <- glm(prediction ~ ns(pred.frame[, 1], df = 8) * 
                      ns(pred.frame[, 2], df = 8), data = pred.frame, family = poisson)
    prediction <- fitted(pred.glm)
  }
  
  max.pred <- max(prediction)
  min.pred <- min(prediction) # new
  message("maximum value = ", round(max.pred, 2), "\n")
  message("minimum value = ", round(min.pred, 2), "\n") # new
  
  if (is.null(z.range)) {
    if (family == "bernoulli") {
      z.range <- c(0, 1)
    } else if (family == "poisson") {
      z.range <- c(0, max.pred * 1.1)
    } else {
      # z.min <- min(data[, y], na.rm = T) ### FEHLER?! Warum z-range aus y?
      # z.max <- max(data[, y], na.rm = T) ### FEHLER?! Warum z-range aus y?
      # z.delta <- z.max - z.min
      # z.range <- c(z.min - (1.1 * z.delta), z.max + (1.1 * z.delta))
      # z.min <- min(data[, y], na.rm = T) 
      # z.max <- max(data[, y], na.rm = T)
      z.delta <- max.pred - min.pred
      z.range <- c(min.pred - (1.1 * z.delta), max.pred + (1.1 * z.delta))
    }
  }
  
  if (have.factor == FALSE) {
    pred.matrix <- matrix(prediction, ncol = 50, nrow = 50)
    if (smooth == "average") {
      pred.matrix.smooth <- pred.matrix
      for (i in 2:49) {
        for (j in 2:49) {
          pred.matrix.smooth[i, j] <- mean(pred.matrix[c((i - 
                                                            1):(i + 1)), c((j - 1):(j + 1))])
        }
      }
      pred.matrix <- pred.matrix.smooth
    }
    
    if (mask) {
      mask.trees <- gbm.object$gbm.call$best.trees
      point.prob <- gbm::predict.gbm(gbm.object[[1]], pred.frame, 
                                     n.trees = mask.trees, type = "response")
      point.prob <- matrix(point.prob, ncol = 50, nrow = 50)
      pred.matrix[point.prob < 0.5] <- 0
    }
    
    if (!perspective) {
      image(x = x.var, y = y.var, z = pred.matrix, zlim = z.range)
    } else {
      persp(x = x.var, y = y.var, z = pred.matrix, zlim = z.range, 
            xlab = x.label, ylab = y.label, zlab = z.label, 
            theta = theta, phi = phi, r = sqrt(10), d = 3, 
            ticktype = ticktype, mgp = c(4, 1, 0), ...)
    }
  }
  
  if (have.factor) {
    factor.list <- names(table(pred.frame[, 1]))
    n <- 1
   
    if (is.null(z.range)) {
      vert.limits <- c(0, max.pred * 1.1)
    } else {
      vert.limits <- z.range
    }
    
    # plot(pred.frame[pred.frame[, 1] == factor.list[1], 2], 
    #      prediction[pred.frame[, 1] == factor.list[1]], type = "l", 
    #      ylim = vert.limits, xlab = y.label, ylab = z.label, 
    #      ...)
    plot(pred.frame[pred.frame[, 1] == factor.list[1], 2], 
         prediction[pred.frame[, 1] == factor.list[1]], type = "l", 
         ylim = vert.limits, xlab = y.label, ylab = z.label)
    
    for (i in 2:length(factor.list)) {
      factor.level <- factor.list[i]
      lines(pred.frame[pred.frame[, 1] == factor.level, 
                       2], prediction[pred.frame[, 1] == factor.level], 
            lty = i)
    }
    
    if (is.null(leg.coords)) {
      x.max <- max(pred.frame[, 2])
      x.min <- min(pred.frame[, 2])
      x.range <- x.max - x.min
      x.pos <- c(x.min + (0.02 * x.range), x.min + (0.3 * 
                                                      x.range))
      y.max <- max(prediction)
      y.min <- min(prediction)
      y.range <- y.max - y.min
      y.pos <- c(y.min + (0.8 * y.range), y.min + (0.95 * 
                                                     y.range))
      legend(x = x.pos, y = y.pos, factor.list, lty = c(1:length(factor.list)), 
             bty = "n")
    } else {
      legend(x = leg.coords[1], y = leg.coords[2], factor.list, 
             lty = c(1:length(factor.list)), bty = "n")
    }
  }
}