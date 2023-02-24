## BRT-tutorial
# https://rspatial.org/raster/sdm/9_sdm_brt.html


# setup -------------------------------------------------------------------
library(gbm)
library(dismo)


# load data ---------------------------------------------------------------
data("Anguilla_train")
head(Anguilla_train)



# identify optimal number of trees (nt) -----------------------------------
angaus.tc5.lr01 <- gbm.step(data            = Anguilla_train,
                            gbm.x           = 3:13,
                            gbm.y           = 2,
                            family          = "bernoulli",
                            tree.complexity = 5,
                            learning.rate   = 0.01,
                            bag.fraction    = 0.5)
summary(angaus.tc5.lr01)


angaus.tc5.lr005 <- gbm.step(data            = Anguilla_train,
                             gbm.x           = 3:13,
                             gbm.y           = 2,
                             family          = "bernoulli",
                             tree.complexity = 5,
                             learning.rate   = 0.005,
                             bag.fraction    = 0.5)
summary(angaus.tc5.lr005)


# plotting ----------------------------------------------------------------
gbm.plot(angaus.tc5.lr005,
         n.plots     = 11,
         plot.layout = c(4, 3),
         write.title = FALSE)

gbm.plot.fits(angaus.tc5.lr005)



# interactions ------------------------------------------------------------
find.int <- gbm.interactions(angaus.tc5.lr005)
find.int$rank.list
find.int$interactions

gbm.perspec(angaus.tc5.lr005,
            x       = 7,
            y       = 1,
            y.range = c(15, 20),
            z.range = c(0, 0.6))
