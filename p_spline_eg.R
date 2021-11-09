
library(mgcv)
library(gratia)
library(tidymv)
library(patchwork)
library(tidyverse)

# NOTE: models done in small sub-sample as memory issues when fitting to whole sample
# will update when can run in whole sample

#### LOAD DATASET #####

dat <- read.csv(
  "dat/alsp_ht_long.csv",
  header=TRUE, 
  sep=",", 
  na.strings=c("","NA"))

dat <- dat %>% select(-age_) %>% 
  rename(age = agey, ht = ht_) %>% 
  select(id, aln, qlet, sex, age, ht) %>% 
  filter(n() > 1 & age>=0.5 & age<=18)%>% 
  drop_na

dat_m <- dat %>% filter(sex=="male") %>% group_by(id) %>% ungroup
dat_m$id_grp <- as.factor(dat_m$id)

dat_f <- dat %>% filter(sex=="female") %>% group_by(id) %>% ungroup
dat_f$id_grp <- as.factor(dat_f$id)

#### P-SPLINE GAMM ####

# MALES

dat_m_small <- dat_m %>% filter(id < 300)

gam_m1 <- gam(
  ht ∼ s(sqrt(age), bs = "ps", k = 7) + 
    s(sqrt(age), id_grp, k = 7, bs = "fs"),
  data = dat_m_small, method = "REML")

ht_m_pred_mean <- predict_gam(
  gam_m1, exclude_terms = "s(sqrt(age),id_grp)")

ht_m_pred <- predict(gam_m1, se.fit = TRUE)

ht_m_pred_dat <- transform(
  dat_m_small, 
  height_fitted = ht_m_pred$fit, 
  height_fitted_se = ht_m_pred$se.fit)

# FEMALES

dat_f_small <- dat_f %>% filter(id < 300)

gam_f1 <- gam(
  ht ∼ s(sqrt(age), bs = "ps", k = 7) + 
    s(sqrt(age), id_grp, k = 7, bs = "fs"),
  data = dat_f_small, method = "REML")

ht_f_pred_mean <- predict_gam(
  gam_f1, exclude_terms = "s(sqrt(age),id_grp)")

ht_f_pred <- predict(gam_f1, se.fit = TRUE)

ht_f_pred_dat <- transform(
  dat_f_small, 
  height_fitted = ht_f_pred$fit, 
  height_fitted_se = ht_f_pred$se.fit)

# PLOT 

(ht_ps_mean_plot_m <- ggplot(
  data = ht_m_pred_mean, aes(
    x = age, y = fit, ymin = fit-(1.96*se.fit), ymax = fit+(1.96*se.fit))) +
    geom_line() + geom_ribbon(alpha=0.3) + theme_classic() + ggtitle(
      "Mean height curve: ALSPAC Males")  + 
    theme(text = element_text(size = 10)))

(indiv_fitted_traj_m <- ggplot(data = ht_m_pred_dat, aes(x = age, y = height_fitted, col = id_grp)) +
    geom_line() + theme_classic() + theme(legend.position = "none") + ggtitle(
      "Individual height curves: ALSPAC Males")  + 
    theme(text = element_text(size = 10)))

(ht_ps_mean_plot_f <- ggplot(
  data = ht_f_pred_mean, aes(
    x = age, y = fit, ymin = fit-(1.96*se.fit), ymax = fit+(1.96*se.fit))) +
    geom_line() + geom_ribbon(alpha=0.3) + theme_classic() + ggtitle(
      "Mean height curve: ALSPAC Females")  + 
    theme(text = element_text(size = 10)))

(indiv_fitted_traj_f <- ggplot(data = ht_f_pred_dat, aes(x = age, y = height_fitted, col = id_grp)) +
    geom_line() + theme_classic() + theme(legend.position = "none") + ggtitle(
      "Individual height curves: ALSPAC Females")  + 
    theme(text = element_text(size = 10)))

graphics.off()

pdf("res/ht_p_spl_plots.pdf",
    width = 7.5, height = 5.5)

((ht_ps_mean_plot_m + indiv_fitted_traj_m) / 
    (ht_ps_mean_plot_f + indiv_fitted_traj_f))

dev.off()

