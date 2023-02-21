# setup -------------------------------------------------------------------
library(tidyverse)
library(betareg)


# logit -------------------------------------------------------------------
library(boot)
tibble(x = seq(0.1, 0.9, 0.1)) %>% 
  mutate(a = logit(x),
         b = inv.logit(a))



# prediction --------------------------------------------------------------
data("FoodExpenditure", package = "betareg")
fe_beta <- betareg(I(food/income) ~ income + persons,
                   data = FoodExpenditure)

summary(fe_beta)

fe2 <- FoodExpenditure %>% 
  mutate(f_i = food / income,
         f_i_pred_aut = predict(fe_beta, .),
         f_i_pred_man = inv.logit(fe_beta$coefficients$mean[1] +
                                    fe_beta$coefficients$mean[2] * income +
                                    fe_beta$coefficients$mean[3] * persons))

rm(fe_beta, fe2, FoodExpenditure)


# link functions ----------------------------------------------------------
tibble(x       = seq(0.01, 0.99, 0.01),
       logit   = log(x / (1 - x)),
       cloglog = log(-log(1 - x)),
       log     = log(x),
       loglog  = -log(-log(x))) %>% 
  pivot_longer(cols = c(logit:loglog),
               names_to = "linkfun") %>% 
  ggplot(aes(x, value, group = linkfun, color = linkfun)) +
  geom_line()

tibble(x       = seq(0.01, 0.99, 0.01),
       logit   = log(x / (1 - x)),
       cloglog = log(-log(1 - x))) %>% 
  pivot_longer(cols = c(logit:cloglog),
               names_to = "linkfun") %>% 
  ggplot(aes(x, value, group = linkfun, color = linkfun)) +
  geom_line()


# mgm int int / 10 --------------------------------------------------------
library(margins)
library(effects)

data_LT_A <- read_rds("data/processed/analysis/analysis_data_transf.rds") %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "A") %>% 
  mutate(mgm_interval_10  = mgm_interval / 10,
         mgm_intensity_10 = mgm_intensity / 10)

m1 <- betareg(sha_i_MP_met_abs ~ stratum + q_site2 * q_reg2 + mgm_type + mgm_interval + mgm_intensity + I(mgm_interval^2) + I(mgm_intensity^2) + mgm_interval * mgm_intensity,
              data = data_LT_A,
              link = "logit")


m2 <- betareg(sha_i_MP_met_abs ~ stratum + q_site2 * q_reg2 + mgm_type + mgm_interval_10 + mgm_intensity_10 + I(mgm_interval_10^2) + I(mgm_intensity_10^2) + mgm_interval_10 * mgm_intensity_10,
              data = data_LT_A,
              link = "logit")

summary(m1)
m1_m1 <- margins(m1)
m1_m2 <- margins(m1, change = "minmax")
e <- allEffects(mod = m1)
plot(allEffects(mod = m1))
summary(m1_m1)
summary(m1_m2)
sm1m2 <- summary(m1_m2) %>% 
  as_tibble()

ggplot(sm1m2, aes(factor, AME)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower,
                    ymax = upper),
                width = 0.25) +
  geom_hline(yintercept = 0) +
  coord_flip()


summary(m2)

plot(allEffects(mod = m2))

rm(data_LT_A, e, m1, m1_m, m2)

# effekte f2 & f3 ---------------------------------------------------------
library(effects)
library(ggeffects)

data_LT_A <- read_rds("data/processed/analysis/analysis_data_transf.rds") %>% 
  filter(simtype == "LT") %>% 
  filter(nat_haz == "A") %>% 
  mutate(mgm_interval_10  = mgm_interval / 10,
         mgm_intensity_10 = mgm_intensity / 10)

f2 <- betareg(sha_i_IP_met_abs ~ stratum + q_site2 * q_reg2 + mgm_type + mgm_interval + mgm_intensity + I(mgm_interval^2) + I(mgm_intensity^2),
              data = data_LT_A,
              link = "logit")
eff_f2 <- allEffects(f2, xlevel = 50)
plot(eff_f2)

eff_f2$mgm_interval$x
eff_f2$`I(mgm_interval^2)`$x
eff_f2$mgm_interval$fit
eff_f2$`I(mgm_interval^2)`$fit

plot(predictorEffect("q_site2", f2))
peff_f2_minterval <- predictorEffect("mgm_interval", f2)

peff_f2_minterval$x
head(peff_f2_minterval$fit)
head(eff_f2$mgm_interval$fit)


ggeff <- ggeffect(f2, c("q_site2", "q_reg2"))
plot(ggeff)

f3 <- betareg(sha_i_IP_met_abs ~ stratum + q_site2 * q_reg2 + mgm_type + mgm_interval * mgm_intensity + I(mgm_interval^2) + I(mgm_intensity^2),
              data = data_LT_A,
              link = "logit")
eff_f3 <- allEffects(f3, xlevel = 50)
plot(eff_f3)

ggeff <- ggeffect(f3, c("mgm_interval", "mgm_intensity"))
plot(ggeff, facets = TRUE)
ggplot(ggeff, aes(x, predicted)) +
  geom_line() +
  facet_wrap(~group)

rm(list = ls())
