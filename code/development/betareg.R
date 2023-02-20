library(tidyverse)
library(boot)
library(betareg)

tibble(x = seq(0.1, 0.9, 0.1)) %>% 
  mutate(a = logit(x),
         b = inv.logit(a))



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
