library(lme4)
library(rpart)
library(tidyverse)

data <- lme4::sleepstudy

# Dataviz
ggplot(data, aes(x = Days, y = Reaction)) +
  geom_point() +
  geom_smooth(color = "pink", size = 0.5, alpha = 0) +
  facet_wrap(~Subject, nrow = 3) +
  theme_classic() +
  scale_x_continuous(labels = scales::pretty_breaks()) +
  labs(x = "Days")

# Tree - can it capture the Subject importance?
tree <- rpart(Reaction ~ Days + Subject, data = data)
tree$variable.importance
rpart.plot::rpart.plot(tree)

# A: yes

m0 <- lm(Reaction ~ Days , data = sleepstudy)
summary(m0)

lmm_days <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy)
lmm <- lmer(Reaction ~ Days + (1 | Subject), sleepstudy)
summary(lmm)
round(confint(lmm, oldNames = FALSE), 2)

anova(lmm, lmm2)


predict_no_re = predict(lmm_days)
predict_lm = predict(m0)


sqrt(sum((data$Reaction - predict_no_re)^2))
sqrt(sum((data$Reaction - predict_lm)^2))

# Bayesian
library(brms)

fit1 <- brm(Reaction ~ Days + (Days | Subject), sleepstudy,
            family = gaussian(),
            prior = c(set_prior("normal(0,5)", class = "b"),
                      set_prior("cauchy(0,2)", class = "sd"),
                      set_prior("lkj(2)", class = "cor")),
            warmup = 1000, iter = 2000, chains = 4,
            control = list(adapt_delta = 0.95))
summary(fit1)

pp <- predict(fit1)
summary(pp)
sqrt(sum((data$Reaction - pp)^2))
